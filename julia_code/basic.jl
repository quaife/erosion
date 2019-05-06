# basic.jl: Basic routines such as datatypes and Stokes solvers.

#--------------- OBJECTS ---------------#
# ThetaLenType: Includes the geometry data and stress of a single body.
mutable struct ThetaLenType
	theta::Vector{Float64}; len::Float64; xsm::Float64; ysm::Float64;
	xx::Vector{Float64}; yy::Vector{Float64}; atau::Vector{Float64}; 
end
# ThLenDenType: Includes the vector of all thlens and the density function.
mutable struct ThLenDenType
	thlenvec::Vector{ThetaLenType}; 
	density::Vector{Float64}; denrot::Vector{Float64};
end
# ParamType: Includes the parameters dt, epsilon, sigma, etc.
mutable struct ParamType
	dt::Float64; epsilon::Float64; sigma::Float64; 
	nouter::Int; ifmm::Int; ibary::Int; maxl::Int; fixarea::Bool; fixpdrop::Bool;
	npts::Int; tfin::Float64; cntout::Int; cput0::Float64;
	geofile::AbstractString; paramsfile::AbstractString
end
# DerivsType: Includes the derivatives of theta, len, xsm, ysm.
mutable struct DerivsType
	mterm::Float64; nterm::Vector{Float64}; 
	xsmdot::Float64; ysmdot::Float64
end
# TargetsType: includes x-y coordinates of target points and u,v,pressure.
mutable struct TargetsType
	xtar::Vector{Float64}; ytar::Vector{Float64};
	utar::Vector{Float64}; vtar::Vector{Float64};
	ptar::Vector{Float64}; vortar::Vector{Float64};
end

#--------------- OBJECT ROUTINES ---------------#
# Create new instances of each type.
function new_thlenden(thlenvec::Vector{ThetaLenType}, density=evec(), denrot=evec())
	return ThLenDenType(thlenvec,density,denrot)
end
function new_thlenvec(nbods::Int)
	return [new_thlen() for nn=1:nbods]
end
function new_thlen()
	return ThetaLenType(evec(),0.,0.,0.,evec(),evec(),evec())
end
function new_dvec(nbods::Int)
	return [new_derivs() for nn=1:nbods]
end
function new_derivs()
	return DerivsType(0.,evec(),0.,0.)
end
function evec()
	return Array{Float64}(undef,0)
end

#--------------- THE MAIN ROUTINE TO GET THE STRESS ---------------#
#= getstress! The main function for calling the necessary Fortran routines.
Computes the smoothed stress atau and saves it in thlenden.thlenvec.atau. =#
function getstress!(thlenden::ThLenDenType, params::ParamType)
	# Compute the density if not loaded already.
	compute_density!(thlenden, params)
	# Compute the stress.
	tau = compute_stress(thlenden,params.nouter,params.ibary,fixpdrop=params.fixpdrop,rotation=false)
	# Smooth atau and save it in each of the thlen variables.
	npts,nbods = getnvals(thlenden.thlenvec)
	for nn = 1:nbods
		n1,n2 = n1n2(npts,nn)
		atau = abs.(tau[n1:n2])
		atau = gaussfilter(atau, params.sigma)
		thlenden.thlenvec[nn].atau = atau[:]
	end
	return
end

#--------------- FORTRAN WRAPPERS ---------------#
#--- THE DENSITY FUNCTION ---#
#= compute_density! Computes the density function and saves in thlenden.
Note: Only computes if density is not already loaded.
Note: It also computes xx and yy along the way and saves in thlenden.thlenvec. =#
function compute_density!(thlenden::ThLenDenType, params::ParamType; rotation::Bool=false)
	if (rotation == false && length(thlenden.density) == 0)
		println("Computing the density function.")
		npts,nbods,xv,yv = getnxy(thlenden)
		thlenden.density = compute_density(xv,yv,npts,nbods,params.nouter,params.ifmm,params.ibary,params.maxl)
	elseif (rotation == true && length(thlenden.denrot) == 0)
		println("Computing the rotated density function.")
		npts,nbods,xv,yv = getnxy(thlenden)
		xrot,yrot = xyrot(xv,yv)
		thlenden.denrot = compute_density(xrot,yrot,npts,nbods,params.nouter,params.ifmm,params.ibary,params.maxl)
	end
	return
end
# compute_density: Fortran wrapper.
function compute_density(xx::Vector{Float64}, yy::Vector{Float64}, 
		npts::Int, nbods::Int, nouter::Int, ifmm::Int, ibary::Int, maxl::Int)
	density = zeros(Float64, 2*npts*nbods + 3*nbods + 2*nouter)
	nits = zeros(Int,1)
	# Call the Fortran routine StokesSolver.
	ccall((:stokessolver_, "libstokes.so"), Nothing, 
		(Ref{Int},Ref{Int},Ref{Int},Ref{Int},Ref{Int},Ref{Int},
		Ref{Float64},Ref{Float64},Ref{Float64},Ref{Int}), 
		npts, nbods, nouter, ifmm, ibary, maxl, xx, yy, density, nits)
	println("The total number of GMRES iterations is ", nits[1],"\n\n")
	return density
end

#--- THE SHEAR STRESS ---#
# compute_stress: Dispatch for ThLenDenType.
function compute_stress(thlenden::ThLenDenType, nouter::Int, ibary::Int; 
		fixpdrop::Bool=false, rotation::Bool=false)
	npts,nbods,xv,yv,density = getnxyden(thlenden,nouter,ibary,fixpdrop,rotation)
	tau = compute_stress(xv,yv,density,npts,nbods,nouter,ibary)
	return tau
end
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

#--- THE PRESSURE ---#
# compute_pressure: Dispatch for ThLenDenType.
function compute_pressure(thlenden::ThLenDenType, nouter::Int, ibary::Int;
		fixpdrop::Bool=false, rotation::Bool=false)
	npts,nbods,xv,yv,density = getnxyden(thlenden,nouter,ibary,fixpdrop,rotation)
	pressure = compute_pressure(xv,yv,density,npts,nbods,nouter,ibary)
	return pressure
end
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
#--- THE QUANTITIES OF INTEREST ---#
# compute_qoi_targets! Dispatch for ThLenDenType and TargetsType. 
function compute_qoi_targets!(thlenden::ThLenDenType, targets::TargetsType, nouter::Int, ibary::Int;
		fixpdrop::Bool=false, rotation::Bool=false)
	npts,nbods,xv,yv,density = getnxyden(thlenden,nouter,ibary,fixpdrop,rotation)
	targets.utar, targets.vtar, targets.ptar, targets.vortar = 
			compute_qoi_targets(xv,yv,density,targets.xtar,targets.ytar,npts,nbods,nouter,ibary)
	return
end
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

#--------------- SMALL ROUTINES ---------------#
# getnxy: For ThLenDenType, get npts, nbods and the x-y coordinates of all the bodies.
function getnxy(thlenden::ThLenDenType)
	npts,nbods = getnvals(thlenden.thlenvec)
	xv,yv = getallxy(thlenden.thlenvec,npts,nbods)
	return npts,nbods,xv,yv
end
# getallxy: Get the x-y coordinates for all of the bodies.
function getallxy(thlenv::Vector{ThetaLenType}, npts::Int, nbods::Int)
	xv,yv = [zeros(Float64,npts*nbods) for ii=1:2]
	for nn = 1:nbods
		getxy!(thlenv[nn])
		n1,n2 = n1n2(npts,nn)
		xv[n1:n2], yv[n1:n2] = thlenv[nn].xx, thlenv[nn].yy
	end
	return xv,yv
end
# getnvals: Calculate npts and nbods.
function getnvals(thlenv::Vector{ThetaLenType})
	nbods = length(thlenv)
	if nbods == 0
		npts = 0
	else
		npts = length(thlenv[1].theta)
	end
	return npts,nbods
end
# Calculate n1 and n2 to divy up the separate bodies.
function n1n2(npts::Integer, nn::Integer)
	n1 = npts*(nn-1)+1
	n2 = npts*nn
	return n1,n2
end
# xyrot: Rotate the x and y coordinates by 90 degrees CCW.
function xyrot(xv::Vector{Float64}, yv::Vector{Float64})
	xrot = -yv
	yrot = xv
	return xrot,yrot
end

#--------------- KEEP PRESSURE DROP FIXED ---------------#
# getnxyden: Get these values depending on fixpdrop and rotation.
function getnxyden(thlenden::ThLenDenType, nouter::Int, ibary::Int, fixpdrop::Bool, rotation::Bool)
	npts,nbods,xv,yv = getnxy(thlenden)
	# Consider fixpdrop.
	rescale = getumax(thlenden, nouter, ibary, fixpdrop)
	# Consider rotation.
	if rotation
		xv,yv = xyrot(xv,yv)
		density = rescale * thlenden.denrot
	else
		density = rescale * thlenden.density
	end
	return npts,nbods,xv,yv,density
end
# getumax: Get umax to rescale the density function.
function getumax(thlenden::ThLenDenType, nouter::Int, ibary::Int, fixpdrop::Bool)
	# NOTE: With u = 1-y^2 and x0 = 2, the pressure drop is pdrop = 8.
	umax = 1.
	if fixpdrop
		pdrop = getpdrop(thlenden, nouter, ibary)[1]
		umax =  10 * 8/pdrop
		println("Fixing pdrop, umax = ", round(umax,sigdigits=3))
	end
	return umax
end
#= getpdrop: Calculate the pressure drop from -x0 to x0. 
Also get the average flux while at it. 
Note: this routine assumes that umax = 1; If different, need to apply rescaling. =#
function getpdrop(thlenden::ThLenDenType, nouter::Int, ibary::Int, x0::Float64 = 2.0; rotation::Bool=false)
	# Set up targets points on two vertical slices.
	nypts = 13
	dy = 2/nypts
	ylocs = collect(-1+0.5*dy: dy: 1-0.5*dy)
	# Target points for plus/minus x0.
	tarp = regulargridtargs([x0],ylocs)
	tarm = regulargridtargs([-x0],ylocs)
	compute_qoi_targets!(thlenden,tarp,nouter,ibary,rotation=rotation)
	compute_qoi_targets!(thlenden,tarm,nouter,ibary,rotation=rotation)
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