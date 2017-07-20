# basic.jl: Basic routines such as datatypes and Stokes solvers.

#--------------- OBJECTS ---------------#
# ParamType: includes the parameters dt, epsilon, sigma, etc.
type ParamType
	dt::Float64; epsilon::Float64; sigma::Float64; nouter::Int; ifmm::Int; fixarea::Int;
end
# ThetaLenType: includes the geometry data for each body and memory terms.
type ThetaLenType
	theta::Vector{Float64}; len::Float64; xsm::Float64; ysm::Float64;
	xx::Vector{Float64}; yy::Vector{Float64};
	atau::Vector{Float64}; mterm::Float64; nterm::Vector{Float64}; 
	xsmdot::Float64; ysmdot::Float64;
end
# ThLenDenType: includes the vector of all thlen's and the density function.
type ThLenDenType
	thlenvec::Vector{ThetaLenType}; density::Vector{Float64}; denrot::Vector{Float64};
end
# TargetsType: includes x-y coordinates of target points and u,v,pressure.
type TargetsType
	xtar::Vector{Float64}; ytar::Vector{Float64};
	utar::Vector{Float64}; vtar::Vector{Float64}; ptar::Vector{Float64};
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
	return ThetaLenType(evec(),0.,0.,0.,evec(),evec(),evec(),0.,evec(),0.,0.)
end
function evec()
	return Array(Float64,0)
end
# Copy all contents from thlen1 to thlen2.
function copy_thlen!(thlen1::ThetaLenType, thlen2::ThetaLenType)
	thlen2.theta = thlen1.theta
	thlen2.len = thlen1.len
	thlen2.xsm = thlen1.xsm
	thlen2.ysm = thlen1.ysm
	thlen2.xx = thlen1.xx
	thlen2.yy = thlen1.yy
	thlen2.atau = thlen1.atau
	thlen2.mterm = thlen1.mterm
	thlen2.nterm = thlen1.nterm
	thlen2.xsmdot = thlen1.xsmdot
	thlen2.ysmdot = thlen1.ysmdot
	return
end

#--------------- MAIN ROUTINES ---------------#
#= getstress! The main function for calling the necessary Fortran routines.
Computes the smoothed stress atau and saves it in thlenden.thlenvec.atau. =#
function getstress!(thlenden::ThLenDenType, params::ParamType)
	# Compute the density (if not loaded already).
	compute_density!(thlenden, params)
	# Compute the stress.
	tau = compute_stress(thlenden, params.nouter)	
	# Smooth atau and save it in each of the thlen variables.
	npts,nbods = getnvals(thlenden.thlenvec)
	for nn = 1:nbods
		n1,n2 = n1n2(npts,nn)
		atau = abs(tau[n1:n2])
		atau = gaussfilter(atau, params.sigma)
		if params.fixarea == 1
			atau = atau - mean(atau)
		end
		thlenden.thlenvec[nn].atau = atau[:]
	end
	return
end
# getrotdensity! Compute the density on the grid rotated by 90 deg CCW.
function getrotdensity!(thlenden::ThLenDenType, params::ParamType)
	npts,nbods,xv,yv = getnxy(thlenden)
	xrot,yrot = xyrot(xv,yv)
	thlenden.denrot = compute_density(xrot,yrot,npts,nbods,params.nouter,params.ifmm)
	return
end

#--------------- FORTRAN WRAPPERS ---------------#
#= computedensity! Computes the density function and saves in thlenden.
Note: Only computes if density is not already loaded.
Note: It also computes xx and yy along the way and saves in thlenden.thlenvec. =#
function compute_density!(thlenden::ThLenDenType, params::ParamType)
	if thlenden.density == []
		npts,nbods,xv,yv = getnxy(thlenden)
		thlenden.density = compute_density(xv,yv,npts,nbods,params.nouter,params.ifmm)
	end
	return
end
# compute_density: Fortran wrapper.
function compute_density(xx::Vector{Float64}, yy::Vector{Float64}, 
		npts::Int, nbods::Int, nouter::Int, ifmm::Int)
	density = zeros(Float64, 2*npts*nbods + 3*nbods + 2*nouter)
	ccall((:stokessolver_, "libstokes.so"), Void, 
		(Ptr{Int},Ptr{Int},Ptr{Int},Ptr{Int},Ptr{Float64},Ptr{Float64},Ptr{Float64}), 
		&npts, &nbods, &nouter, &ifmm, xx, yy, density)
	return density
end

# computestress: Dispatch for ThLenDenType.
function compute_stress(thlenden::ThLenDenType, nouter::Int)
	npts,nbods,xv,yv = getnxy(thlenden)
	if nbods == 0
		tau = evec()
	else
		tau = compute_stress(xv,yv,thlenden.density,npts,nbods,nouter)
	end
	return tau
end
# compute_stress: Fortran wrapper.
function compute_stress(xx::Vector{Float64}, yy::Vector{Float64}, density::Vector{Float64}, 
		npts::Int, nbods::Int, nouter::Int)
	tau = zeros(Float64, npts*nbods)
	ccall((:computeshearstress_, "libstokes.so"), Void,
		(Ptr{Int},Ptr{Int},Ptr{Int},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64}),
		&npts, &nbods, &nouter, xx, yy, density, tau)
	return tau
end

# computepressure: Dispatch for ThLenDenType.
function compute_pressure(thlenden::ThLenDenType, nouter::Int)
	npts,nbods,xv,yv = getnxy(thlenden)
	if nbods ==0
		pressure = evec()
	else
		pressure = compute_pressure(xv,yv,thlenden.density,npts,nbods,nouter)
	end
	return pressure
end
# compute_pressure: Fortran wrapper.
function compute_pressure(xx::Vector{Float64}, yy::Vector{Float64}, density::Vector{Float64}, 
		npts::Int, nbods::Int, nouter::Int)
	pressure = zeros(Float64, npts*nbods)
	ccall((:computepressure_, "libstokes.so"), Void,
		(Ptr{Int},Ptr{Int},Ptr{Int},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64}),
		&npts, &nbods, &nouter, xx, yy, density, pressure)
	return pressure
end

# compute_velpress_targets: Dispatch for ThLenDenType and TargetsType.
function compute_velpress_targets!(thlenden::ThLenDenType, targets::TargetsType, nouter::Int)
	npts,nbods,xv,yv = getnxy(thlenden)
	targets.utar,targets.vtar,targets.ptar = compute_velpress_targets(xv,yv,
		thlenden.density,targets.xtar,targets.ytar,npts,nbods,nouter)
	return
end
# compute_velpress_targets: Fortran wrapper.
function compute_velpress_targets(xx::Vector{Float64}, yy::Vector{Float64},
		density::Vector{Float64}, xtar::Vector{Float64}, ytar::Vector{Float64},
		npts::Int, nbods::Int, nouter::Int)
	ntargets = endof(xtar)
	utar,vtar,ptar = [zeros(Float64,ntargets) for ii=1:3]
	ccall((:computevelocitypressuretargets_, "libstokes.so"), Void,
		(Ptr{Int},Ptr{Int},Ptr{Int},Ptr{Float64},Ptr{Float64},Ptr{Float64},
		Ptr{Int}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
		&npts, &nbods, &nouter, xx, yy, density, 
		&ntargets, xtar, ytar, utar, vtar, ptar)
	return utar,vtar,ptar
end

#--------------- OTHER ---------------#
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
	if nbods ==0
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
end
