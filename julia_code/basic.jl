# basic.jl: Basic routines such as datatypes and Stokes solvers.

#--------------- OBJECTS ---------------#
# ThetaLenType: Includes the geometry data and stress of a single body.
mutable struct ThetaLenType
	theta::Vector{Float64}; len::Float64; xsm::Float64; ysm::Float64; matau::Float64
end
# ThLenDenType: Includes the vector of all thlens and the density function.
mutable struct ThLenDenType
	thlenvec::Vector{ThetaLenType}; tt::Float64;
	density::Vector{Float64}; denrot::Vector{Float64}; 
end
new_thlenden(thlenvec::Vector{ThetaLenType}) = ThLenDenType(thlenvec, 0, [],[])
# TargetsType: includes x-y coordinates of target points and u,v,pressure.
mutable struct TargetsType
	xtar::Vector{Float64}; ytar::Vector{Float64};
	utar::Vector{Float64}; vtar::Vector{Float64};
	ptar::Vector{Float64}; vortar::Vector{Float64};
end

#--------------- THE MAIN ROUTINE TO GET THE STRESS ---------------#
#= getstress! The main function for calling the necessary Fortran routines.
Computes the smoothed stress atau and saves it in thlenden.thlenvec.atau. =#
function getstress!(thlenden::ThLenDenType, params::ParamSet)
	# Compute the density and stress.
	dummy, den_time = @timed( compute_density!(thlenden, params) )
	stress, str_time = @timed( compute_stress(thlenden,params,fixpdrop=params.fixpdrop) )
	println("\nTime taken to compute density = ", round(den_time,sigdigits=3), "sec.")
	println("Time taken to compute stress = ", round(str_time,sigdigits=3), "sec.")
	# Smooth the stress and also save the mean of smoothed stress atau.
	for bod = 1:length(thlenden.thlenvec)
		stress[:,bod] = gaussfilter( abs.(stress[:,bod]), params.sigma)
		thlenden.thlenvec[bod].matau = mean(stress[:,bod])
	end
	return stress
end

#--------------- FORTRAN WRAPPERS ---------------#
#--- THE DENSITY FUNCTION ---#
# compute_density: Fortran wrapper.
function compute_density(xx::Vector{Float64}, yy::Vector{Float64}, nbods::Int, params::ParamSet)
	@unpack npts, nouter, ifmm, ibary, ibc, maxl = params
	density = zeros(Float64, 2*npts*nbods + 3*nbods + 2*nouter)
	nits = zeros(Int,1)
	# Call the Fortran routine StokesSolver.
	ccall((:stokessolver_, "libstokes.so"), Nothing, 
		(Ref{Int},Ref{Int},Ref{Int},Ref{Int},Ref{Int},Ref{Int},Ref{Int},
		Ref{Float64},Ref{Float64},Ref{Float64},Ref{Int}), 
		npts, nbods, nouter, ifmm, ibary, ibc, maxl, xx, yy, density, nits)
	println("The total number of GMRES iterations is ", nits[1],"\n\n")
	return density
end
# compute_density! Computes the density function and saves in thlenden.
function compute_density!(thlenden::ThLenDenType, params::ParamSet; rotation::Bool=false)
	density = rotation ? thlenden.denrot : thlenden.density
	if length(density) > 0; return; end
	println("Computing the density function with rotation = ", rotation)
	nbods,xv,yv = getnxy(thlenden)
	xv,yv = rotation ? xyrot(xv,yv) : (xv,yv)
	density = compute_density(xv,yv,nbods,params)
	rotation ? thlenden.denrot = density : thlenden.density = density 
end

#--- THE SHEAR STRESS ---#
# compute_stress: Fortran wrapper.
function compute_stress(xx::Vector{Float64}, yy::Vector{Float64}, 
		density::Vector{Float64}, npts::Int, nbods::Int, nouter::Int, ibary::Int)
	tau = zeros(Float64, npts*nbods)
	if nbods > 0
		ccall((:computeshearstress_, "libstokes.so"), Nothing,
			(Ref{Int},Ref{Int},Ref{Int},
			Ref{Float64},Ref{Float64},Ref{Float64},Ref{Int},Ref{Float64}),
			npts, nbods, nouter, xx, yy, density, ibary, tau)
	end
	return tau
end
# compute_stress: Dispatch for ThLenDenType.
function compute_stress(thlenden::ThLenDenType, params::ParamSet; 
		fixpdrop::Bool=false, rotation::Bool=false)
	@unpack npts, nouter, ibary = params
	nbods,xv,yv,density = getnxyden(thlenden,params,fixpdrop,rotation)
	tau = compute_stress(xv,yv,density,npts,nbods,nouter,ibary)
	return reshape(tau,npts,nbods)
end

#--- THE PRESSURE ---#
# compute_pressure: Fortran wrapper.
function compute_pressure(xx::Vector{Float64}, yy::Vector{Float64}, 
		density::Vector{Float64}, npts::Int, nbods::Int, nouter::Int, ibary::Int)
	pressure = zeros(Float64, npts*nbods)
	if nbods > 0
		ccall((:computepressure_, "libstokes.so"), Nothing,
			(Ref{Int},Ref{Int},Ref{Int},
			Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64}),
			npts, nbods, nouter, xx, yy, density, pressure)
	end
	return pressure
end
# compute_pressure: Dispatch for ThLenDenType.
function compute_pressure(thlenden::ThLenDenType, params::ParamSet;
		fixpdrop::Bool=false, rotation::Bool=false)
	@unpack npts, nouter, ibary = params
	nbods,xv,yv,density = getnxyden(thlenden,params,fixpdrop,rotation)
	pressure = compute_pressure(xv,yv,density,npts,nbods,nouter,ibary)
	return pressure
end

#--- THE QUANTITIES OF INTEREST ---#
# compute_qoi_targets: Fortran wrapper.
function compute_qoi_targets(xx::Vector{Float64}, yy::Vector{Float64},
		density::Vector{Float64}, xtar::Vector{Float64}, ytar::Vector{Float64},
		npts::Int, nbods::Int, nouter::Int, ibary::Int)
	ntargets = length(xtar)
	utar,vtar,ptar,vortar = [zeros(Float64,ntargets) for ii=1:4]
	ccall((:computeqoitargets_, "libstokes.so"), Nothing,
		(Ref{Int}, Ref{Int}, Ref{Int}, Ref{Int}, Ref{Float64}, Ref{Float64}, Ref{Float64},
		Ref{Int}, Ref{Float64}, Ref{Float64}, 
		Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
		npts, nbods, nouter, ibary, xx, yy, density, 
		ntargets, xtar, ytar, utar, vtar, ptar, vortar)
	return utar,vtar,ptar,vortar
end
# compute_qoi_targets! Dispatch for ThLenDenType and TargetsType. 
function compute_qoi_targets!(thlenden::ThLenDenType, targets::TargetsType, params::ParamSet;
		fixpdrop::Bool=false, rotation::Bool=false)
	@unpack npts, nouter, ibary = params
	nbods,xv,yv,density = getnxyden(thlenden,params,fixpdrop,rotation)
	targets.utar, targets.vtar, targets.ptar, targets.vortar = 
			compute_qoi_targets(xv,yv,density,targets.xtar,targets.ytar,npts,nbods,nouter,ibary)
	return
end


#--------------- SMALL ROUTINES ---------------#
# xyrot: Rotate the x and y coordinates by 90 degrees CCW.
xyrot(xv::Vector{Float64}, yv::Vector{Float64}) = -yv, xv
# getnxy: Get npts, nbods, and the x-y coordinates of all the bodies.
function getnxy(thlenden::ThLenDenType)
	nbods = length(thlenden.thlenvec)
	xv,yv = [Array{Float64}(undef,0) for ii=1:2]
	for nn = 1:nbods
		xx,yy = getxy(thlenden.thlenvec[nn])
		append!(xv,xx); append!(yv,yy)
	end
	return nbods,xv,yv
end
# getnxyden: Get these values depending on fixpdrop and rotation.
function getnxyden(thlenden::ThLenDenType, params::ParamSet, fixpdrop::Bool, rotation::Bool)
	nbods,xv,yv = getnxy(thlenden)
	# Consider fixpdrop.
	rescale = getumax(thlenden, params, fixpdrop)
	# Consider rotation.
	if rotation
		xv,yv = xyrot(xv,yv)
		density = rescale * thlenden.denrot
	else
		density = rescale * thlenden.density
	end
	return nbods,xv,yv,density
end

#--------------- KEEP PRESSURE DROP FIXED ---------------#
# getumax: Get umax to rescale the density function.
function getumax(thlenden::ThLenDenType, params::ParamSet, fixpdrop::Bool, x0::Float64 = 2.0)
	# NOTE: With u = 1-y^2 and x0 = 2, the pressure drop is pdrop = 8.
	# This code should work for pipe flow or uniform flow: ibc = 0 or 1.
	umax = 1.
	if fixpdrop
		pdrop = getpdrop(thlenden, params)[1]
		umax = 4*x0/(pdrop + 1e-4)
		println("Fixing pdrop, umax = ", round(umax,sigdigits=3))
	end
	return umax
end
#= getpdrop: Calculate the pressure drop from -x0 to x0. 
Also get the average flux while at it. 
Note: this routine assumes that umax = 1; If different, need to apply rescaling. =#
function getpdrop(thlenden::ThLenDenType, params::ParamSet, x0::Float64 = 2.0; rotation::Bool=false)
	# Set up targets points on two vertical slices.
	nypts = 13
	dy = 2/nypts
	ylocs = collect(-1+0.5*dy: dy: 1-0.5*dy)
	# Target points for plus/minus x0.
	tarp = regulargridtargs([x0],ylocs)
	tarm = regulargridtargs([-x0],ylocs)
	compute_qoi_targets!(thlenden,tarp,params,rotation=rotation)
	compute_qoi_targets!(thlenden,tarm,params,rotation=rotation)
	# Compute the pressure drop.
	pplus = mean(tarp.ptar)
	pminus = mean(tarm.ptar)
	pdrop = pminus-pplus
	# Compute the flux.
	qplus = mean(tarp.utar)
	qminus = mean(tarm.utar)
	qavg = 0.5*(qplus+qminus)
	#= The discharge should be exactly the same at any location x.
	So check that it is the same at x0 and -x0. =#
	qreldiff = (qplus-qminus)/qavg
	if qreldiff > 1e-3
		@warn string("The flux does not match at x0 and -x0: qreldiff = ", round(qreldiff,sigdigits=3))
	end
	return pdrop,qavg
end
# regulargrid: Set up target points on a regular grid; return x and y.
function regulargrid(xlocs::Vector{Float64}, ylocs::Vector{Float64})
	nx = length(xlocs)
	ny = length(ylocs)
	ntargs = nx*ny
	xtar = zeros(Float64,ntargs)
	ytar = zeros(Float64,ntargs)
	for nn=1:nx
		n1 = ny*(nn-1)+1
		n2 = ny*nn
		xtar[n1:n2] .= xlocs[nn]
		ytar[n1:n2] .= ylocs
	end
	return xtar,ytar
end
# regulargridtargs: Set up target points on a regular grid; return targets.
function regulargridtargs(xlocs::Vector{Float64}, ylocs::Vector{Float64})
	xtar,ytar = regulargrid(xlocs,ylocs)
	targets = TargetsType([], [], [], [], [], [])
	targets.xtar = xtar
	targets.ytar = ytar
	return targets
end