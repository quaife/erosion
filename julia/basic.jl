# basic.jl: Basic routines such as datatypes and Stokes solvers.

#################### Object data types ####################
# ParamType: includes the parameters dt, epsilon, sigma, etc.
type ParamType
	dt::Float64; epsilon::Float64; sigma::Float64; ifmm::Int; fixarea::Int;
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
	thlenvec::Vector{ThetaLenType}; density::Vector{Float64};
end
#################### Object routines #####################
# Create a new ThetaLenType that has all zeros.
function new_thlen()
	return ThetaLenType([], 0., 0., 0., [], [], [], 0., [], 0., 0.)
end
# evec: Create an empty Vector{Float64}
function evec()
	return Array(Float64,0)
end
##########################################################

#################### Call Fortran routines ####################
#= getstress! Computes the smoothed stress atau and saves it in thlenden.thlenvec.atau.
Note: It computes the density function only if not done already. =#
function getstress!(thlenden::ThLenDenType, params::ParamType)
	# Compute the density (if not done already).
	getdensity!(thlenden, params.ifmm)
	# Will use the thlenvec component.
	thlenv = thlenden.thlenvec
	npts,nbods,ntot = npnb(thlenv)
	# Compute the stress.
	tau = getstress(xv,yv,density,npts,nbods,ntot)
	# Smooth atau and save it in each of the thlen variables.
	for nn = 1:nbods
		n1,n2 = n1n2(npts,nn)
		atau = abs(tau[n1:n2])
		atau = gaussfilter(atau, params.sigma)
		if params.fixarea == 1; atau = atau - mean(atau); end
		thlenv[nn].atau = atau
	end
	return
end
# getstress: Wrapper for Fortran routine 'computeShearStress' to compute the shear stress.
function getstress(xx::Vector{Float64}, yy::Vector{Float64}, density::Vector{Float64},
		npts::Int, nbods::Int, ntot::Int)
	tau = zeros(Float64, ntot)
	ccall((:computeShearStress_, "libstokes.so"), Void,
		(Ptr{Int},Ptr{Int},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64}),
		&npts, &nbods, xx, yy, density, tau)
	return tau
end
#= getdensity! Computes the density function (if not done already) and saves in thlenden.
Note: It also computes xx and yy along the way and saves in thlenden.thlenvec. =#
function getdensity!(thlenden::ThLenDenType, ifmm::Int)
	if thlenden.density == []
		thlenv = thlenden.thlenvec
		npts,nbods,ntot = npnb(thlenv)
		xv,yv = getallxy(thlenv,npts,nbods,ntot)
		thlenv.density = getdensity(xv,yv,npts,nbods,ntot,ifmm)
	end
	return
end
# getdensity: Wrapper for Fortran routine 'stokesSolver' to get the density function.
function getdensity(xx::Vector{Float64}, yy::Vector{Float64}, 
		npts::Int, nbods::Int, ntot::Int, ifmm::Int)
	density = zeros(Float64, ntot)
	ccall((:stokessolver_, "libstokes.so"), Void, 
		(Ptr{Int},Ptr{Int},Ptr{Int},Ptr{Float64},Ptr{Float64},Ptr{Float64}), 
		&npts, &nbods, &ifmm, xx, yy, density)
	return density
end
# getallxy: Get the x-y coordinates for all of the bodies.
function getallxy(thlenv::Vector{ThetaLenType},npts::Int,nbods::Int,ntot::Int)
	xv,yv = [zeros(Float64,ntot) for ii=1:2]
	for nn = 1:nbods
		getxy!(thlenv[nn])
		n1,n2 = n1n2(npts,nn)
		xv[n1:n2], yv[n1:n2] = thlenv[nn].xx, thlenv[nn].yy
	end
	return xv,yv
end

# Small routines 
# npnb: Calculate npts and nbods
function npnb(thlenvec::Vector{ThetaLenType})
	nbods = endof(thlenv)
	npts = endof(thlenv[1].theta)
	ntot = npts*nbods
	return npts,nbods,ntot
end
# Calculate n1 and n2 to divy up the separate bodies.
function n1n2(npts::Integer,nn::Integer)
	n1 = npts*(nn-1)+1
	n2 = npts*nn
	return n1,n2
end
##################################################
