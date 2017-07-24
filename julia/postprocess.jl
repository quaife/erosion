# postporcess.jl
# Post-processing routines.

# postprocess: Use the saved data to compute stuff.
function postprocess(foldername::AbstractString)
	# Define the data folder.
	datafolder = string("../datafiles/",foldername,"/")
	# Read the params data file.
	paramsfile = string(datafolder,"params.dat")
	paramvec = readvec(paramsfile)
	nouter = paramvec[2]
	ntimes = paramvec[end]
	# Read the data at each time step.
	for cnt=0:ntimes
		# Get the file name at each time.
		cntstr = lpad(cnt,4,0)
		geomfile = string(datafolder,"geom",cntstr,".dat")
		densityfile = string(datafolder,"density",cntstr,".dat")
		# Extract thlenvec, density, and denrot.
		tt,thlenvec = read_geom_file(geomfile)
		density,denrot = read_density_file(densityfile)
		# Create variables.
		npts,nbods = getnvals(thlenvec)
		thlenden = new_thlenden(thlenvec,density,denrot)

		#--------------------------------------#
		# Compute the resistivity (1/permeability) of the matrix.
		rbods = resistivity(thlenden, nouter, 2.0, false)
		rbodsrot = resistivity(thlenden, nouter, 2.0, true)
		# Compute the drag on each body.
		dragx, dragy = drag(thlenden, nouter, false)
		dragxrot, dragyrot = drag(thlenden, nouter, true)
		# Save the data to a file.
		resdragfile = string(datafolder,"resdrag",cntstr,".dat")
		lab1 = string("# Data on resistivity and drag: ")
		lab2 = string("# nbods, resistivity, rotated resistivity, 
			total dragx and dragy, rotated dragx and dragy.")
		resdragdata = [lab1; lab2; nbods; rbods; rbodsrot; 
			dragx; dragy; dragxrot; dragyrot]
		writedata(resdragdata, resdragfile)

		#--------------------------------------#
		# Compute velocity and pressure at a set of target points.
		xlocs = [-2.5, -2.0, -1.5, 1.5, 2.0, 2.5]
		ylocs = collect(-0.8: 0.2: 0.8)
		targets = setuptargets(xlocs,ylocs)
		compute_velpress_targets!(thlenden,targets,nouter)
		# Save the output to a data file.
		targfile = string(datafolder,"targs",cntstr,".dat")
		label = string("# Data at grid of target points: x, y, u, v, pressure.")
		targdata = [label; targets.xtar; targets.ytar; 
			targets.utar; targets.vtar; targets.ptar]
		writedata(targdata, targfile)
	end
	return
end

# resistivity: Compute the resistivity/permeability of the porous matrix.
function resistivity(thlenden::ThLenDenType, nouter::Int, x0::Float64, rotation::Bool=false)
	# Set up targets points on a y-grid for midpoint rule.
	nypts = 13
	dy = 2./nypts
	ylocs = collect(-1+0.5*dy: dy: 1-0.5*dy)
	# Target points for plus/minus x0.
	tarp = setuptargets([x0],ylocs)
	tarm = setuptargets([-x0],ylocs)
	# Compute the velocities and pressures on each set of target points.
	# Either using the original porous matrix or the rotated one.
	if rotation==true
		compute_velpressrot_targets!(thlenden,tarp,nouter)
		compute_velpressrot_targets!(thlenden,tarm,nouter)
	else
		compute_velpress_targets!(thlenden,tarp,nouter)
		compute_velpress_targets!(thlenden,tarm,nouter)
	end
	# Compute the cross-sectional average pressure and discharge
	pplus = mean(tarp.ptar)
	pminus = mean(tarm.ptar)
	qplus = mean(tarp.utar)
	qminus = mean(tarm.utar)
	qavg = 0.5*(qplus+qminus)
	#= The discharge should be exactly the same at any location x.
	So check that it is the same at x0 and -x0. =#
	qreldiff = (qplus-qminus)/qavg
	assert(qreldiff < 1e-6)
	# Calculate the total resistivity
	rtot = (pminus - pplus)/(2*x0*qavg)
	# Calculate the resisitvity due only to the bodies.
	rbods = x0*(rtot - 3)
	# For testing.
	#println("At x0 = ", x0, " the total resistivity is: ", signif(rtot,3))
	#println("At x0 = ", x0, " the matrix resistivity is: ", signif(rbods,3))
	return rbods
end
# drag: Compute the total drag on all of the bodies combined.
function drag(thlenden::ThLenDenType, nouter::Int, rotation::Bool=false)
	# Get the shear stress and pressure on the set of bodies.
	# Note: the stress is not smoothed and absolute value is not taken.
	if rotation==true
		tauvec = compute_stressrot(thlenden,nouter)
		pressvec = compute_pressrot(thlenden,nouter)
		
		### CHECK THIS
		thlenvec = thlenden.thlenvec + 0.5*pi
		###

	else
		tauvec = compute_stress(thlenden,nouter)
		pressvec = compute_pressure(thlenden,nouter)
		thlenvec = thlenden.thlenvec
	end
	npts,nbods = getnvals(thlenvec)
	dragx = 0.
	dragy = 0.
	for nn=1:nbods
		# Get the pressure and stress and body nn.
		n1,n2 = n1n2(npts,nn)
		press = pressvec[n1:n2]
		tau = tauvec[n1:n2]
		# Get the tangent/normal vectors and arc length increment.
		sx,sy,nx,ny = getns(thlenvec[nn].theta)
		ds = thlenvec[nn].len / npts
		# Compute the drag force.
		# Note: I believe both should be plus signs due to the conventions of s and n.
		dragx += sum(press.*nx + tau.*sx)*ds
		dragy += sum(press.*ny + tau.*sy)*ds
	end
	return dragx, dragy
end
# setuptargets: Set up the target points.
function setuptargets(xlocs::Vector{Float64}, ylocs::Vector{Float64})
	nx = endof(xlocs)
	ny = endof(ylocs)
	ntargs = nx*ny
	xtar = ones(Float64,ntargs)
	ytar = ones(Float64,ntargs)
	for nn=1:nx
		n1,n2 = n1n2(ny,nn)
		xtar[n1:n2] = xlocs[nn]
		ytar[n1:n2] = ylocs
	end
	targets = TargetsType(evec(), evec(), evec(), evec(), evec())
	targets.xtar = xtar
	targets.ytar = ytar
	return targets
end
