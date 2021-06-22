# MAIN GOAL: Postprocess the jld2 file to compute new quantities.


using LinearAlgebra
include("main.jl")

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
		reldiff > 1e-3 ? @warn(string("Relative error in area = ", round(reldiff,sigdigits=2))) : 0.
		areavec[nn] = area
	end
	return areavec
end
# resistivity: Compute the resistivity/permeability of the porous matrix.
function resistivity(thlenden::ThLenDenType, params::ParamType, x0::Float64=2.0; rotation::Bool=false)
	# Retrieve the pressure drop and flux (assuming umax = 1)
	pdrop,qavg = getpdrop(thlenden,params.nouter,params.ibary,x0,rotation=rotation)	
	# Calculate the total resistivity
	resist = pdrop/(2*x0*qavg)
	# If pipe flow (ibc = 0) then remove the contribution from the walls.
	if params.ibc == 0; resist = x0*(resist - 3); end
	return resist
end
# drag: Compute the total drag on all of the bodies combined.
function drag(thlenden::ThLenDenType, params::ParamType; rotation::Bool=false)
	# Get the shear stress and pressure on the set of bodies.
	# Note 1: the stress is not smoothed and absolute value is not taken.
	# Note 2: these are the values with umax = 1.
	tauvec = compute_stress(thlenden,params.nouter,params.ibary,fixpdrop=false,rotation=rotation)
	pressvec = compute_pressure(thlenden,params.nouter,params.ibary,fixpdrop=false,rotation=rotation)
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
		atauvec[n1:n2] = gaussfilter(abs.(tau), params.sigma)
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
	umax = getumax(thlenden, params.nouter, params.ibary,params.fixpdrop)
	return pdragx, pdragy, vdragx, vdragy, umax, tauvec, atauvec
end

#----------- ROUTINES FOR TARGET POINT CALCULATIONS -----------#
# regbodtargs: Set up target points on a regular and body fitted grid.
function regbodtargs(thlenv::Vector{ThetaLenType})
	# Regular grid.
	hh = 0.05
	xlocs = collect(-1-2*hh: hh: 1+2*hh)
	ylocs = collect(-1+0.5*hh: hh: 1-0.5*hh)	
	xreg,yreg = regulargrid(xlocs,ylocs)
	# Body-fitted grid.
	spacevec = 0.01*collect(1:2:3)
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
	nlayers = length(spacevec)
	xtar,ytar = Array{Float64}(undef,0), Array{Float64}(undef,0)
	for nn = 1:nbods
		thlen = thlenv[nn]
		xx,yy = thlen.xx[ind], thlen.yy[ind]
		nx,ny = getns(thlen.theta)[3:4]
		nx,ny = nx[ind], ny[ind]
		append!(xtar, vec(xx*ones(1,nlayers) - nx*transpose(spacevec)))
		append!(ytar, vec(yy*ones(1,nlayers) - ny*transpose(spacevec)))
		# Remove the points that lie outside the computational domain.
		badind = findall( abs.(ytar) .> 0.999 )
		deleteat!(xtar,badind)
		deleteat!(ytar,badind)
	end
	return xtar,ytar
end




# getns: Get the normal and tangent directions.
# Convention: CCW parameterization and inward pointing normal.
function getns(theta::Vector{Float64}, rotation::Bool=false)
	# CCW tangent vector.
	if rotation == false
		sx, sy = cos.(theta), sin.(theta)
	else
		sx, sy = -sin.(theta), cos.(theta)
	end
	# Inward pointing normal vector.
	nx, ny = -sy, sx
	return sx,sy,nx,ny
end






#----------- MAIN ROUTINES -----------#

#----------- STARTUP AND ASSISTING ROUTINES -----------#
# startpostprocess
function startpostprocess(foldername::AbstractString)
	# Define the data folder and files.
	datafolder = string("../output_data/", foldername, "/")
	paramsfile = string(datafolder, "aparams")
	pinfofile = string(datafolder, "apinfo.txt")
	# Get extra information from apinfo.
	pinfovec = readvec(pinfofile)
	ntimes = Int(pinfovec[2])
	# Get the params object.
	params = getparams(paramsfile)
	return datafolder, ntimes, params
end


# pp1: Postprocess the fast stuff: area and resistivity.
function pp1(datafile::AbstractString)
	println("Beginning pp1 on ", datafile)
	
	#datafolder,ntimes,params = startpostprocess(datafile)
	## Need ntime params


	# Read the data at each time step.
	for cnt=0:ntimes
		print("pp1 step ", cnt, " of ", ntimes, ": ")
		##thlenden = GET
		npts, nbods = getnvals(thlenden.thlenvec)
		# Compute the area of each body.
		areavec = getareas(thlenden)
		if nbods == 0 areavec = [0] end
		print("area completed; ")
		# Compute the resistivity (1/permeability) of the matrix.
		rbods = resistivity(thlenden,params)
		rbodsrot = resistivity(thlenden,params,rotation=true)
		println("resistivity completed; ")
	end


		# Save the data to a file.
		## Write areavec to file
		add_data(datafile, areavec)
		# Save the resistivity data to the same JLD file.

		iofile = jldopen(file, "r+")
		write(iofile, "y", y)
		close(iofile)



	println("Finished pp1 on ", foldername)
	return
end




# pp2: Postprocess the slower stuff: drag and stress.
function pp2(foldername::AbstractString)
	println("\n\nBeginning pp2 on ", foldername)
	datafolder,ntimes,params = startpostprocess(foldername)
	# Read the data at each time step.
	for cnt=0:ntimes
		println("\npp2 beggining step ", cnt, " of ", ntimes, ".")
		thlenden, cntstr = get_thlenden(datafolder,cnt)
		npts,nbods = getnvals(thlenden.thlenvec)
		#--------------------------------------#
		# Compute the total pressure and viscous drag on the collection of bodies.
		pdrx,pdry,vdrx,vdry,umax,tauv,atauv = drag(thlenden,params)
		pdrxr,pdryr,vdrxr,vdryr,umaxr,tauvr,atauvr = drag(thlenden,params,rotation=true)
		println("Finished the drag computation.")
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
	end
	println("Finished pp2 on ", foldername, "\n")
end


# pp3: Postprocess the slowest stuff: quantities of interest at the target points.
function pp3(foldername::AbstractString)
	println("Beginning pp3 on ", foldername)
	datafolder,ntimes,params = startpostprocess(foldername)
	# Read the data at each time step.
	for cnt=0:ntimes
		print("pp3 step ", cnt, " of ", ntimes, "; ")
		thlenden, cntstr = get_thlenden(datafolder,cnt)
		npts,nbods = getnvals(thlenden.thlenvec)
		#--------------------------------------#
		# Compute velocity, pressure, vorticity at a set of target points, with umax set to 1.
		targets = regbodtargs(thlenden.thlenvec)
		compute_qoi_targets!(thlenden,targets,params.nouter,params.ibary,fixpdrop=false)
		# Save the output to a data file.
		targfile = string(datafolder,"targs",cntstr,".dat")
		label = string("# Data at grid of target points: x, y, u, v, pressure, vorticity.")
		targdata = [label; targets.xtar; targets.ytar; 
			targets.utar; targets.vtar; targets.ptar; targets.vortar]
		writedata(targdata, targfile, digs=5)
		println("step completed")
	end
	println("Finished pp3 on ", foldername, "\n")
end



#postprocess: Run all postprocess routines pp1-3.
function postprocess(datafile::AbstractString)
	println("\n\n%------------------------------------------------------%")
	t1 = time()
	println("Beginning postprocessing ", datafile, "\n")
	pp1(datafile)
#	pp2(datafile)
#	pp3(datafile)
	println("Finished postprocessing ", datafile)
	pptime = round((time()-t1)/60., sigdigits=2)
	println("Time taken: ", pptime, " minutes.")
	println("%------------------------------------------------------%\n\n")
end