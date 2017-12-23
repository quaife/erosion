# postporcess.jl
# Post-processing routines.

# postprocess: Use the saved data to compute stuff.
function postprocess(foldername::AbstractString)
	# Define the data folder and files.
	datafolder = string("../datafiles/",foldername,"/")
	paramsfile = string(datafolder,"aparams.in")
	pinfofile = string(datafolder,"apinfo.out")
	# Get extra information from apinfo.
	pinfovec = readvec(pinfofile)
	npts = Int(pinfovec[1])
	ntimes = Int(pinfovec[3])
	# Get the params object.
	params = getparams(paramsfile,npts)
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
		# Compute the area of each body.
		areavec = getareas(thlenden)
		# Save the data to a file.
		areasfile = string(datafolder,"areas",cntstr,".dat")
		label = string("# Area of each individual body.")
		areadata = [label; areavec]
		writedata(areadata, areasfile)

		#--------------------------------------#
		# Compute the resistivity (1/permeability) of the matrix.
		rbods = resistivity(thlenden, nouter, 2.0)
		rbodsrot = resistivity(thlenden, nouter, 2.0, rotation=true)
		# Compute the total pressure and viscous drag on the collection of bodies.
		pdragx, pdragy, vdragx, vdragy, tauvec = drag(thlenden, params)
		pdragxr, pdragyr, vdragxr, vdragyr, tauvecr = drag(thlenden, params, rotation=true)
		# Save the data to a file.
		resdragfile = string(datafolder,"resdrag",cntstr,".dat")
		lab1 = string("# Data on resistivity and drag: ")
		lab2 = string("# nbods, resistivity, rotated resistivity, pdragx, pdragy, vdragx, vdragy, and all rotated values.")
		resdragdata = [lab1; lab2; nbods; rbods; rbodsrot; 
			pdragx; pdragy; vdragx; vdragy; pdragxr; pdragyr; vdragxr; vdragyr]
		writedata(resdragdata, resdragfile)

		#--------------------------------------#
		# Compute velocity, pressure, vorticity at a set of target points.
		targets = regbodtargs(thlenvec)
		compute_qoi_targets!(thlenden,targets,nouter, fixpdrop = params.fixpdrop)
		# Save the output to a data file.
		targfile = string(datafolder,"targs",cntstr,".dat")
		label = string("# Data at grid of target points: x, y, u, v, pressure, vorticity.")
		targdata = [label; targets.xtar; targets.ytar; 
			targets.utar; targets.vtar; targets.ptar; targets.vortar]
		writedata(targdata, targfile)

		#--------------------------------------#
		# Save the stress on each body.
		# atauvec has absolute value and smoothing applied; 
		# tauvec is raw stress, with a + or - sign and no smoothing.
		getstress!(thlenden,params)
		stressfile = string(datafolder,"stress",cntstr,".dat")
		label = string("# Smoothed atau, Raw atau ")
		atauvec = zeros(Float64, npts*nbods)
		for nn=1:nbods
			n1,n2 = n1n2(npts,nbods)
			atauvec[n1:n2] = thlenden.thlenvec[nn].atau
		end
		writedata([label; atauvec; tauvec], stressfile)
		# Print progress.
		println("Finished step ", cnt, " of ", ntimes, ".")
	end
	return
end

# getareas: Compute the area of each body.
function getareas(thlenden::ThLenDenType)
	npts,nbods = getnvals(thlenden.thlenvec)
	areavec = zeros(Float64,nbods)
	for nn = 1:nbods
		thlen = thlenden.thlenvec[nn]
		xx,yy = thlen.xx,thlen.yy
		len = thlen.len
		sx,sy,nx,ny = getns(thlen.theta)
		# Compute area in two ways to estimate error.
		areax = dot(xx,nx)*len
		areay = dot(yy,ny)*len
		area = 0.5*(areax+areay)
		reldiff = abs(areax-areay)/area
		reldiff > 1e-3? warn("Relative error in area = ", signif(reldiff,2)) : 0.
		areavec[nn] = area
	end
	return areavec
end
# resistivity: Compute the resistivity/permeability of the porous matrix.
function resistivity(thlenden::ThLenDenType, nouter::Int, x0::Float64=2.0; 
		rotation::Bool=false)
	pdrop,qavg = getpdrop(thlenden,nouter,x0, rotation=rotation)
	# Calculate the total resistivity
	rtot = pdrop/(2*x0*qavg)
	# Calculate the resisitvity due only to the bodies.
	rbods = x0*(rtot - 3)
	# For testing.
	#println("At x0 = ", x0, " the total resistivity is: ", signif(rtot,3))
	#println("At x0 = ", x0, " the matrix resistivity is: ", signif(rbods,3))
	return rbods
end
# drag: Compute the total drag on all of the bodies combined.
function drag(thlenden::ThLenDenType, params::ParamType; rotation::Bool=false)
	# Get the shear stress and pressure on the set of bodies.
	# Note: the stress is not smoothed and absolute value is not taken.
	tauvec = compute_stress(thlenden,params.nouter,fixpdrop=params.fixpdrop,rotation=rotation)
	pressvec = compute_pressure(thlenden,params.nouter,fixpdrop=params.fixpdrop,rotation=rotation)
	thlenvec = thlenden.thlenvec
	npts,nbods = getnvals(thlenvec)
	pdragx,pdragy,vdragx,vdragy = 0.,0.,0.,0.
	for nn=1:nbods
		# Get the pressure and stress on body nn.
		n1,n2 = n1n2(npts,nn)
		press = pressvec[n1:n2]
		tau = tauvec[n1:n2]
		# Get the tangent/normal vectors and arc length increment.
		sx,sy,nx,ny = getns(thlenvec[nn].theta,rotation)
		ds = thlenvec[nn].len / npts
		# Compute the drag force.
		# Note: I believe both should be plus signs due to the conventions of s and n.
		# Compute the pressure and viscous drag separately.
		pdragx += sum(press.*nx)*ds
		pdragy += sum(press.*ny)*ds
		vdragx += sum(tau.*sx)*ds
		vdragy += sum(tau.*sy)*ds
	end
	return pdragx, pdragy, vdragx, vdragy, tauvec
end

#----------- TARGET POINTS -----------#
# regbodtargs: Set up target points on a regular and body fitted grid.
function regbodtargs(thlenv::Vector{ThetaLenType})
	# Regular grid.
	hh = 0.10
	xlocs = collect(-2.: hh: 2.)
	ylocs = collect(-1+0.5*hh: hh: 1-0.5*hh)	
	xreg,yreg = regulargrid(xlocs,ylocs)
	# Body-fitted grid.
	spacevec = 0.02*collect(1:2:5)
	nptslayer = 32
	xbod,ybod = bodyfitgrid(thlenv, spacevec, nptslayer)
	# Combine the regular and body fitted grid into a single set of points.
	targets = TargetsType(evec(), evec(), evec(), evec(), evec(), evec())
	targets.xtar = [xreg; xbod]
	targets.ytar = [yreg; ybod]
	return targets
end
# bodyfitgrid: Set up target points on a body fitted grid.
function bodyfitgrid(thlenv::Vector{ThetaLenType}, 
		spacevec::Vector{Float64}, nptslayer::Int)
	npts,nbods = getnvals(thlenv)
	# Use nptslayer in each layer.
	ind0 = div(npts,2*nptslayer)
	ind0 = max(ind0,1)
	ind = ind0:2*ind0:npts
	nlayers = endof(spacevec)
	xtar,ytar = Array(Float64,0), Array(Float64,0)
	for nn = 1:nbods
		thlen = thlenv[nn]
		xx,yy = thlen.xx[ind], thlen.yy[ind]
		nx,ny = getns(thlen.theta)[3:4]
		nx,ny = nx[ind],ny[ind]
		append!(xtar, vec(xx*ones(1,nlayers) - nx*transpose(spacevec)))
		append!(ytar, vec(yy*ones(1,nlayers) - ny*transpose(spacevec)))
	end
	return xtar,ytar
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
# getns: Get the normal and tangent directions.
# Convention: CCW parameterization and inward pointing normal.
function getns(theta::Vector{Float64}, rotation::Bool=false)
	# CCW tangent vector.
	if rotation == false
		sx, sy = cos(theta), sin(theta)
	else
		sx, sy = -sin(theta), cos(theta)
	end
	# Inward pointing normal vector.
	nx, ny = -sy, sx
	return sx,sy,nx,ny
end
