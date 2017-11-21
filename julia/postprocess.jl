# postporcess.jl
# Post-processing routines.

# postprocess: Use the saved data to compute stuff.
function postprocess(foldername::AbstractString)
	# Define the data folder.
	datafolder = string("../datafiles/",foldername,"/")
	# Read the params data file.
	paramsfile = string(datafolder,"params.dat")
	paramvec = readvec(paramsfile)
	npts = paramvec[10]
	ntimes = paramvec[12]
	params = getparams(paramvec[1:9],npts)
	nouter = params.nouter
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
		rbods = resistivity(thlenden, nouter, 2.0)
		rbodsrot = resistivity(thlenden, nouter, 2.0; rotation=true)
		# Compute the drag on each body.
		dragx, dragy = drag(thlenden, nouter)
		dragxrot, dragyrot = drag(thlenden, nouter; rotation=true)
		# Save the data to a file.
		resdragfile = string(datafolder,"resdrag",cntstr,".dat")
		lab1 = string("# Data on resistivity and drag: ")
		lab2 = string("# nbods, resistivity, rotated resistivity, total dragx and dragy, rotated dragx and dragy.")
		resdragdata = [lab1; lab2; nbods; rbods; rbodsrot; 
			dragx; dragy; dragxrot; dragyrot]
		writedata(resdragdata, resdragfile)

		#--------------------------------------#
		# Compute velocity and pressure at a set of target points.
		xlocs = collect(-2.00: 0.05: 2.00)
		ylocs = collect(-0.95: 0.05: 0.95)
		targets = regbodtargs(xlocs,ylocs,thlenvec)
		compute_qoi_targets!(thlenden,targets,nouter)
		# Save the output to a data file.
		targfile = string(datafolder,"targs",cntstr,".dat")
		label = string("# Data at grid of target points: x, y, u, v, pressure, vorticity.")
		targdata = [label; targets.xtar; targets.ytar; 
			targets.utar; targets.vtar; targets.ptar; targets.vortar]
		writedata(targdata, targfile)
		# Print progress.
		println("Finished step ", cnt, " of ", ntimes, ".")

		#--------------------------------------#
		# Save the stress on each body.
		getstress!(thlenden,params)
		stressfile = string(datafolder,"stress",cntstr,".dat")
		label = string("# Data ")
		atauvec = zeros(Float64, npts*nbods)
		for nn=1:nbods
			n1,n2 = n1n2(npts,nbods)
			atauvec[n1:n2] = thlenden.thlenvec[nn].atau
		end
		writedata([label; atauvec], stressfile)
	end
	return
end

# resistivity: Compute the resistivity/permeability of the porous matrix.
function resistivity(thlenden::ThLenDenType, nouter::Int, x0::Float64; rotation::Bool=false)
	# Set up targets points on a y-grid for midpoint rule.
	nypts = 13
	dy = 2./nypts
	ylocs = collect(-1+0.5*dy: dy: 1-0.5*dy)
	# Target points for plus/minus x0.
	tarp = regulargridtargs([x0],ylocs)
	tarm = regulargridtargs([-x0],ylocs)
	# Compute the velocities and pressures on each set of target points.
	# Either using the original porous matrix or the rotated one.
	if rotation==true
		compute_qoirot_targets!(thlenden,tarp,nouter)
		compute_qoirot_targets!(thlenden,tarm,nouter)
	else
		compute_qoi_targets!(thlenden,tarp,nouter)
		compute_qoi_targets!(thlenden,tarm,nouter)
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
	#assert(qreldiff < 1e-6)
	if qreldiff > 1e-6
		warn("The flux does not match at x0 and -x0: qreldiff = ", qreldiff)
	end
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
function drag(thlenden::ThLenDenType, nouter::Int; rotation::Bool=false)
	# Get the shear stress and pressure on the set of bodies.
	# Note: the stress is not smoothed and absolute value is not taken.
	if rotation==true
		tauvec = compute_stressrot(thlenden,nouter)
		pressvec = compute_pressrot(thlenden,nouter)
	else
		tauvec = compute_stress(thlenden,nouter)
		pressvec = compute_pressure(thlenden,nouter)
	end
	thlenvec = thlenden.thlenvec
	npts,nbods = getnvals(thlenvec)
	dragx = 0.
	dragy = 0.
	for nn=1:nbods
		# Get the pressure and stress on body nn.
		n1,n2 = n1n2(npts,nn)
		press = pressvec[n1:n2]
		tau = tauvec[n1:n2]
		# Get the tangent/normal vectors and arc length increment.
		if rotation==true
			sx,sy,nx,ny = getnsrot(thlenvec[nn].theta)
		else
			sx,sy,nx,ny = getns(thlenvec[nn].theta)
		end
		ds = thlenvec[nn].len / npts
		# Compute the drag force.
		# Note: I believe both should be plus signs due to the conventions of s and n.
		dragx += sum(press.*nx + tau.*sx)*ds
		dragy += sum(press.*ny + tau.*sy)*ds
	end
	return dragx, dragy
end

# regbodtargs: Set up target points on a regular and body fitted grid.
function regbodtargs(xlocs::Vector{Float64}, ylocs::Vector{Float64}, thlenv::Vector{ThetaLenType})
	xreg,yreg = regulargrid(xlocs,ylocs)
	xbod,ybod = bodyfitgrid(thlenv)
	targets = TargetsType(evec(), evec(), evec(), evec(), evec(), evec())
	targets.xtar = [xreg; xbod]
	targets.ytar = [yreg; ybod]
	return targets
end
# regulargridtargs: Set up target points on a regular grid; return targets.
function regulargridtargs(xlocs::Vector{Float64}, ylocs::Vector{Float64})
	xtar,ytar = regulargrid(xlocs,ylocs)
	targets = TargetsType(evec(), evec(), evec(), evec(), evec(), evec())
	targets.xtar = xtar
	targets.ytar = ytar
	return targets
end
# regulargrid: Set up target points on a regular grid; return x and y.
function regulargrid(xlocs::Vector{Float64}, ylocs::Vector{Float64})
	nx = endof(xlocs)
	ny = endof(ylocs)
	ntargs = nx*ny
	xtar = zeros(Float64,ntargs)
	ytar = zeros(Float64,ntargs)
	for nn=1:nx
		n1,n2 = n1n2(ny,nn)
		xtar[n1:n2] = xlocs[nn]
		ytar[n1:n2] = ylocs
	end
	return xtar,ytar
end
# bodyfitgrid: Set up target points on a body fitted grid.
function bodyfitgrid(thlenv::Vector{ThetaLenType})
	npts,nbods = getnvals(thlenv)
	xtar,ytar = Array(Float64,0), Array(Float64,0)
	spacevec = transpose(collect(1:2:9))
	mm = endof(spacevec)

	# ROUGH FOR NOW
	nstep = round(Int,npts/64)
	nstep = max(nstep,1)
	ind = nstep:2*nstep:npts

	for nn = 1:nbods
		thlen = thlenv[nn]
		xx,yy = thlen.xx[ind], thlen.yy[ind]
		sx,sy,nx,ny = getns(thlen.theta)
		nxx,nyy = nx[ind], ny[ind]

		#ds = thlen.len/npts
		append!(xtar, vec(xx*ones(1,mm) - 0.02*nxx*spacevec))
		append!(ytar, vec(yy*ones(1,mm) - 0.02*nyy*spacevec))
	end
	return xtar,ytar
end


# getns: Get the normal and tangent directions.
# Convention: CCW parameterization and inward pointing normal.
function getns(theta::Vector{Float64})
	# CCW tangent vector.
	sx, sy = cos(theta), sin(theta)
	# Inward pointing normal vector.
	nx, ny = -sy, sx
	return sx,sy,nx,ny
end
# getns: Get the normal and tangent directions on the rotated grid.
function getnsrot(theta::Vector{Float64})
	# CCW tangent vector.
	sx, sy = -sin(theta), cos(theta)
	# Inward pointing normal vector.
	nx, ny = -sy, sx
	return sx,sy,nx,ny
end
