# basic.jl: Basic routines such as datatypes and Stokes solvers.

#--------------- OBJECTS ---------------#
# ThetaLenType: Includes the geometry data and stress of a single body.
type ThetaLenType
	theta::Vector{Float64}; len::Float64; xsm::Float64; ysm::Float64;
	xx::Vector{Float64}; yy::Vector{Float64}; atau::Vector{Float64}; 
end
# ThLenDenType: Includes the vector of all thlens and the density function.
type ThLenDenType
	thlenvec::Vector{ThetaLenType}; 
	density::Vector{Float64}; denrot::Vector{Float64};
end
# ParamType: Includes the parameters dt, epsilon, sigma, etc.
type ParamType
	dt::Float64; epsilon::Float64; sigma::Float64; 
	nouter::Int; ifmm::Int; fixarea::Bool; fixpdrop::Bool;
	npts::Int; tfin::Float64; cntout::Int; cput0::Float64
end
# DerivsType: Includes the derivatives of theta, len, xsm, ysm
type DerivsType
	mterm::Float64; nterm::Vector{Float64}; 
	xsmdot::Float64; ysmdot::Float64
end
# TargetsType: includes x-y coordinates of target points and u,v,pressure.
type TargetsType
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
	return Array(Float64,0)
end

#--------------- THE MAIN ROUTINE TO GET THE STRESS ---------------#
#= getstress! The main function for calling the necessary Fortran routines.
Computes the smoothed stress atau and saves it in thlenden.thlenvec.atau. =#
function getstress!(thlenden::ThLenDenType, params::ParamType)
	# Compute the density if not loaded already.
	compute_density!(thlenden, params)
	# Compute the stress.
	tau = compute_stress(thlenden, params.nouter, 
		fixpdrop = params.fixpdrop, rotation = false)
	# Smooth atau and save it in each of the thlen variables.
	npts,nbods = getnvals(thlenden.thlenvec)
	for nn = 1:nbods
		n1,n2 = n1n2(npts,nn)
		atau = abs(tau[n1:n2])
		atau = gaussfilter(atau, params.sigma)
		if params.fixarea
			atau = atau - mean(atau)
		end
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
	if (rotation == false & endof(thlenden.density) == 0)
		println("Computing the density function.")
		npts,nbods,xv,yv = getnxy(thlenden)
		thlenden.density = compute_density(xv,yv,npts,nbods,params.nouter,params.ifmm)
	elseif (rotation == true & endof(thlenden.density) == 0)
		println("Computing the rotated density function.")
		npts,nbods,xv,yv = getnxy(thlenden)
		xrot,yrot = xyrot(xv,yv)
		thlenden.denrot = compute_density(xrot,yrot,npts,nbods,params.nouter,params.ifmm)
	end
	return
end
# compute_density: Fortran wrapper.
function compute_density(xx::Vector{Float64}, yy::Vector{Float64}, 
		npts::Int, nbods::Int, nouter::Int, ifmm::Int)
	density = zeros(Float64, 2*npts*nbods + 3*nbods + 2*nouter)
	# Call the Fortran routine StokesSolver.
	ccall((:stokessolver_, "libstokes.so"), Void, 
		(Ptr{Int},Ptr{Int},Ptr{Int},Ptr{Int},Ptr{Float64},Ptr{Float64},Ptr{Float64}), 
		&npts, &nbods, &nouter, &ifmm, xx, yy, density)
	return density
end

#--- THE SHEAR STRESS ---#
# compute_stress: Dispatch for ThLenDenType.
function compute_stress(thlenden::ThLenDenType, nouter::Int; 
		fixpdrop::Bool=false, rotation::Bool=false)
	npts,nbods,xv,yv,density = getnxyden(thlenden,nouter,fixpdrop,rotation)
	tau = compute_stress(xv,yv,density,npts,nbods,nouter)
	return tau
end
# compute_stress: Fortran wrapper.
function compute_stress(xx::Vector{Float64}, yy::Vector{Float64}, 
		density::Vector{Float64}, npts::Int, nbods::Int, nouter::Int)
	tau = zeros(Float64, npts*nbods)
	ccall((:computeshearstress_, "libstokes.so"), Void,
		(Ptr{Int},Ptr{Int},Ptr{Int},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64}),
		&npts, &nbods, &nouter, xx, yy, density, tau)
	return tau
end

#--- THE PRESSURE ---#
# compute_pressure: Dispatch for ThLenDenType.
function compute_pressure(thlenden::ThLenDenType, nouter::Int;
		fixpdrop::Bool=false, rotation::Bool=false)
	npts,nbods,xv,yv,density = getnxyden(thlenden,nouter,fixpdrop,rotation)
	pressure = compute_pressure(xv,yv,density,npts,nbods,nouter)
	return pressure
end
# compute_pressure: Fortran wrapper.
function compute_pressure(xx::Vector{Float64}, yy::Vector{Float64}, 
		density::Vector{Float64}, npts::Int, nbods::Int, nouter::Int)
	pressure = zeros(Float64, npts*nbods)
	ccall((:computepressure_, "libstokes.so"), Void,
		(Ptr{Int},Ptr{Int},Ptr{Int},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64}),
		&npts, &nbods, &nouter, xx, yy, density, pressure)
	return pressure
end
#--- THE QUANTITIES OF INTEREST ---#
# compute_qoi_targets! Dispatch for ThLenDenType and TargetsType. 
function compute_qoi_targets!(thlenden::ThLenDenType, targets::TargetsType, nouter::Int;
		fixpdrop::Bool=false, rotation::Bool=false)
	npts,nbods,xv,yv,density = getnxyden(thlenden,nouter,fixpdrop,rotation)
	targets.utar, targets.vtar, targets.ptar, targets.vortar = 
			compute_qoi_targets(xv,yv,density,targets.xtar,targets.ytar,npts,nbods,nouter)
	return
end
# compute_qoi_targets: Fortran wrapper.
function compute_qoi_targets(xx::Vector{Float64}, yy::Vector{Float64},
		density::Vector{Float64}, xtar::Vector{Float64}, ytar::Vector{Float64},
		npts::Int, nbods::Int, nouter::Int)
	ntargets = endof(xtar)
	utar,vtar,ptar,vortar = [zeros(Float64,ntargets) for ii=1:4]
	ccall((:computeqoitargets_, "libstokes.so"), Void,
		(Ptr{Int}, Ptr{Int}, Ptr{Int}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
		Ptr{Int}, Ptr{Float64}, Ptr{Float64}, 
		Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
		&npts, &nbods, &nouter, xx, yy, density, 
		&ntargets, xtar, ytar, utar, vtar, ptar, vortar)
	return utar,vtar,ptar,vortar
end

#--------------- SMALL ROUTINES ---------------#
# getnxyden: Get these values depending on fixpdrop and rotation.
function getnxyden(thlenden::ThLenDenType, nouter::Int, 
		fixpdrop::Bool, rotation::Bool)
	npts,nbods,xv,yv = getnxy(thlenden)
	# Consider fixpdrop.
	rescale = 1.
	if fixpdrop
		pdrop = getpdrop(thlenden, nouter)[1]
		rescale = 8./pdrop
		# NOTE: With u = 1-y^2 and x0 = 2, the pressure drop is pdrop = 8.
	end
	# Consider rotation.
	if rotation
		xv,yv = xyrot(xv,yv)
		density = rescale * thlenden.denrot
	else
		density = rescale * thlenden.density
	end
	return npts,nbods,xv,yv,density
end
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
	nbods = endof(thlenv)
	if nbods == 0
		npts = 0
	else
		npts = endof(thlenv[1].theta)
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
#= getpdrop: Calculate the pressure drop from -x0 to x0. 
Also get the average flux while at it. =#
function getpdrop(thlenden::ThLenDenType, nouter::Int, 
		x0::Float64 = 2.0, rotation::Bool=false)
	# Set up targets points on two vertical slices.
	nypts = 13
	dy = 2./nypts
	ylocs = collect(-1+0.5*dy: dy: 1-0.5*dy)
	# Target points for plus/minus x0.
	tarp = regulargridtargs([x0],ylocs)
	tarm = regulargridtargs([-x0],ylocs)
	compute_qoi_targets!(thlenden,tarp,nouter, rotation=rotation)
	compute_qoi_targets!(thlenden,tarm,nouter, rotation=rotation)
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
		warn("The flux does not match at x0 and -x0: qreldiff = ", qreldiff)
	end
	return pdrop,qavg
end
