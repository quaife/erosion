# postporcess.jl
# Post-processing routines.

#----------- MAIN ROUTINES -----------#
#postprocess: Run all postprocess routines pp1-3.
function postprocess(foldername::AbstractString)
	pp1(foldername)
	pp2(foldername)
	pp3(foldername)
end

# pp1: Postprocess the fast stuff: area and resistivity.
function pp1(foldername::AbstractString)
	println("Beggining pp1.")
	datafolder,ntimes,params = startpostprocess(foldername)
	# Read the data at each time step.
	for cnt=0:ntimes
		println("Beggining step ", cnt, " of ", ntimes, ".")
		thlenden, cntstr = get_thlenden(datafolder,cnt)
		npts,nbods = getnvals(thlenden.thlenvec)
		#--------------------------------------#
		# Compute the area of each body.
		println("Beginning the area computation.")
		areavec = getareas(thlenden)
		# Save the data to a file.
		areasfile = string(datafolder,"areas",cntstr,".dat")
		label = string("# Area of each individual body.")
		areadata = [label; areavec]
		writedata(areadata, areasfile)
		println("Finished the area computation.\n")
		#--------------------------------------#
		# Compute the resistivity (1/permeability) of the matrix.
		println("Beginning the resistivity computation.")
		rbods = resistivity(thlenden,params.nouter,2.0)
		rbodsrot = resistivity(thlenden,params.nouter,2.0,rotation=true)
		# Save the resistivity data to a file.
		resfile = string(datafolder,"resistivity",cntstr,".dat")
		label = string("# Data on resistivity: nbods, resistivity, rotated resistivity.")
		resdata = [label; nbods; rbods; rbodsrot;]
		writedata(resdata, resfile)
		# Print progress.
		println("Finished the resistivity computation.\n")
		println("Finished step ", cnt, " of ", ntimes, ".\n\n")
	end
	println("Finished pp1.")
	return
end
# pp2: Postprocess the slower stuff: drag and stress.
function pp2(foldername::AbstractString)
	println("Beggining pp2.")
	datafolder,ntimes,params = startpostprocess(foldername)
	# Read the data at each time step.
	for cnt=0:ntimes
		println("Beggining step ", cnt, " of ", ntimes, ".")
		thlenden, cntstr = get_thlenden(datafolder,cnt)
		npts,nbods = getnvals(thlenden.thlenvec)
		#--------------------------------------#
		# Compute the total pressure and viscous drag on the collection of bodies.
		println("Beginning the drag computation.")
		pdrx,pdry,vdrx,vdry,umax,tauv,atauv = drag(thlenden,params)
		pdrxr,pdryr,vdrxr,vdryr,umaxr,tauvr,atauvr = drag(thlenden,params,rotation=true)
		println("Finished the drag computation.\n")
		# Save the data to a file.
		dragfile = string(datafolder,"drag",cntstr,".dat")
		lab1 = string("# Data on drag: ")
		lab2 = string("# nbods, pdragx, pdragy, vdragx, vdragy, umax and all rotated values.")
		dragdata = [lab1; lab2; nbods; 
			pdrx; pdry; vdrx; vdry; umax; pdrxr; pdryr; vdrxr; vdryr]
		writedata(dragdata, dragfile)
		#--------------------------------------#
		# Save the stress on each body.
		# tauv is raw stress, with nontrivial sign and no smoothing.
		# atauv has absolute value and smoothing applied; 
		stressfile = string(datafolder,"stress",cntstr,".dat")
		label = string("# Smoothed atau, Raw atau ")
		writedata([label; atauv; tauv], stressfile)
		println("Finished step ", cnt, " of ", ntimes, ".\n\n")
	end
	println("Finished pp2.")
end
# pp3: Postprocess the slowest stuff: quantities of interest at the target points.
function pp3(foldername::AbstractString)
	println("Beggining pp3.")
	datafolder,ntimes,params = startpostprocess(foldername)
	# Read the data at each time step.
	for cnt=0:ntimes
		println("Beggining step ", cnt, " of ", ntimes, ".")
		thlenden, cntstr = get_thlenden(datafolder,cnt)
		npts,nbods = getnvals(thlenden.thlenvec)
		#--------------------------------------#
		# Compute velocity, pressure, vorticity at a set of target points.
		targets = regbodtargs(thlenden.thlenvec)
		compute_qoi_targets!(thlenden,targets,params.nouter,fixpdrop=params.fixpdrop)
		# Save the output to a data file.
		targfile = string(datafolder,"targs",cntstr,".dat")
		label = string("# Data at grid of target points: x, y, u, v, pressure, vorticity.")
		targdata = [label; targets.xtar; targets.ytar; 
			targets.utar; targets.vtar; targets.ptar; targets.vortar]
		writedata(targdata, targfile)
		println("Finished step ", cnt, " of ", ntimes, ".\n\n")
	end
	println("Finished pp3.")
end


#----------- STARTUP AND ASSISTING ROUTINES -----------#
# startpostprocess
function startpostprocess(foldername::AbstractString)
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
	return datafolder, ntimes, params
end
# get_thlenden
function get_thlenden(datafolder::AbstractString,cnt::Int)
	# Get the file name at each time.
	cntstr = lpad(cnt,4,0)
	geomfile = string(datafolder,"geom",cntstr,".dat")
	densityfile = string(datafolder,"density",cntstr,".dat")
	# Extract thlenvec, density, and denrot.
	tt,thlenvec = read_geom_file(geomfile)
	density,denrot = read_density_file(densityfile)
	# Create variable thlenden
	thlenden = new_thlenden(thlenvec,density,denrot)
	return thlenden, cntstr
end

#----------- ROUTINES FOR AREA, RESISTIVITY, DRAG, ETC. -----------#
# getareas: Compute the area of each body.
function getareas(thlenden::ThLenDenType)
	npts,nbods = getnvals(thlenden.thlenvec)
	areavec = zeros(Float64,nbods)
	for nn = 1:nbods
		thlen = thlenden.thlenvec[nn]
		xx,yy = thlen.xx,thlen.yy
		sx,sy,nx,ny = getns(thlen.theta)
		ds = thlen.len / npts
		# Compute area in two ways to estimate error.
		areax = -dot(xx,nx)*ds
		areay = -dot(yy,ny)*ds
		area = 0.5*(areax+areay)
		reldiff = abs(areax-areay)/area
		reldiff > 1e-3? warn("Relative error in area = ", signif(reldiff,2)) : 0.
		areavec[nn] = area
	end
	return areavec
end
# resistivity: Compute the resistivity/permeability of the porous matrix.
function resistivity(thlenden::ThLenDenType, nouter::Int, x0::Float64=2.0; rotation::Bool=false)
	# Retrieve the pressure drop and flux (assuming umax = 1)
	pdrop,qavg = getpdrop(thlenden,nouter,x0,rotation=rotation)
	# Calculate the total resistivity
	rtot = pdrop/(2*x0*qavg)
	# Calculate the resisitvity due only to the bodies.
	rbods = x0*(rtot - 3)
	return rbods
end

# drag: Compute the total drag on all of the bodies combined.
function drag(thlenden::ThLenDenType, params::ParamType; rotation::Bool=false)
	# Get the shear stress and pressure on the set of bodies.
	# Note 1: the stress is not smoothed and absolute value is not taken.
	# Note 2: these are the values with umax = 1.
	tauvec = compute_stress(thlenden,params.nouter,fixpdrop=false,rotation=rotation)
	pressvec = compute_pressure(thlenden,params.nouter,fixpdrop=false,rotation=rotation)
	thlenvec = thlenden.thlenvec
	npts,nbods = getnvals(thlenvec)
	pdragx,pdragy,vdragx,vdragy = 0.,0.,0.,0.
	atauvec = zeros(Float64,npts*nbods)
	for nn=1:nbods
		# Get the pressure and stress on body nn.
		n1,n2 = n1n2(npts,nn)
		press = pressvec[n1:n2]
		tau = tauvec[n1:n2]
		# Go ahead and compute the absolute smoothed stress.
		atauvec[n1:n2] = gaussfilter(abs(tau), params.sigma)
		# Get the tangent/normal vectors and arc length increment.
		sx,sy,nx,ny = getns(thlenvec[nn].theta,rotation)
		ds = thlenvec[nn].len / npts
		# Compute the pressure and viscous drag separately.
		# Note: I believe both should be plus signs due to the conventions of s and n.
		pdragx += dot(press,nx)*ds
		pdragy += dot(press,ny)*ds
		vdragx += dot(tau,sx)*ds
		vdragy += dot(tau,sy)*ds
	end
	umax = getumax(thlenden, params.nouter, params.fixpdrop)
	return pdragx, pdragy, vdragx, vdragy, umax, tauvec, atauvec
end

#----------- ROUTINES FOR TARGET POINT CALCULATIONS -----------#
# regbodtargs: Set up target points on a regular and body fitted grid.
function regbodtargs(thlenv::Vector{ThetaLenType})
	# Regular grid.
	hh = 0.05
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
