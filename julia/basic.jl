# basic.jl: Basic routines such as datatypes and Stokes solvers.

#################### Object data types ####################
# ParamType includes the parameters dt, epsilon, and beta.
type ParamType
	dt::Float64; epsilon::Float64; sigma::Float64; beta::Real;
end
# ThetaLenType includes all of the data that for a curve.
type ThetaLenType
	theta::Vector{Float64}; len::Float64; xsm::Float64; ysm::Float64;
	xx::Vector{Float64}; yy::Vector{Float64};
	atau::Vector{Float64}; mterm::Float64; nterm::Vector{Float64}; 
	xsmdot::Float64; ysmdot::Float64;
end
##################################################

#################### Includes ####################
using Winston
include("spectral.jl")
include("thetalen.jl")
include("misc.jl")
include("main.jl")
##################################################

#################### Object routines ####################
# Create a new ThetaLenType that has all zeros.
function new_thlen()
	return ThetaLenType([], 0., 0., 0., [], [], [], 0., [], 0., 0.)
end
# Copy the relevant contents from thlen1 to thlen2.
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
# evec: Create an empty Vector{Float64}
function evec()
	return Array(Float64,0)
end
##################################################

#################### Stokes solvers ####################
# All routines work for multiple bodies.
# stokes: Julia wrapper to call the Fortran stokessolver
function stokes(npts::Integer, nbods::Integer, xx::Vector{Float64}, yy::Vector{Float64},
		ntargs::Integer, xtar::Vector{Float64}, ytar::Vector{Float64})
	ntot = npts*nbods
	tau = zeros(Float64, ntot)
	utar = zeros(Float64, ntargs)
	vtar = zeros(Float64, ntargs)
	ptar = zeros(Float64, ntargs)
	ccall((:stokessolver_, "libstokes.so"),
		Void, (Ptr{Int}, Ptr{Int}, Ptr{Int}, 
		Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
		Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}), 
		&npts, &nbods, &ntargs, xx, yy, xtar, ytar, utar, vtar, ptar, tau)
	return tau, utar, vtar, ptar
end
#= stokes!: Dispatch for vector of ThetaLenType
Calculates atau = abs(tau) and smooths it with a Gaussian filter;
then loads each atau in thlenv. =#
function stokes!(thlenv::Vector{ThetaLenType}, sigma::Float64, 
			ntargs::Integer=0, xtar::Vector{Float64}=evec(), ytar::Vector{Float64}=evec())
	nbods = endof(thlenv)
	npts = endof(thlenv[1].theta)
	ntot = nbods*npts
	xv,yv = [zeros(Float64,ntot) for ii=1:2]
	# Put all of the xy values in a single vector.
	for nn = 1:nbods
		getxy!(thlenv[nn])
		n1,n2 = n1n2(npts,nn)
		xv[n1:n2], yv[n1:n2] = thlenv[nn].xx, thlenv[nn].yy
	end
	# Call the stokessolver.
	tau,utar,vtar,ptar = stokes(npts,nbods,xv,yv,ntargs,xtar,ytar)
	# Smooth atau and save it in each of the thlen variables.
	for nn = 1:nbods
		n1,n2 = n1n2(npts,nn)
		atau = abs(tau[n1:n2])
		atau = gaussfilter(atau,sigma)
		thlenv[nn].atau = atau
	end
	return utar,vtar,ptar
end
# Calculate n1 and n2.
function n1n2(npts::Integer,nn::Integer)
	n1 = npts*(nn-1)+1
	n2 = npts*nn
	return n1,n2
end
##################################################
