# basic.jl: Basic stuff.

#################### Object data types ####################
# ParamType includes the parameters dt, epsilon, and beta.
type ParamType
	dt::Float64; epsilon::Float64; beta::Real;
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
##################################################

#################### Object functions ####################
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
##################################################

#################### Stokes functions ####################
# All routines work for multiple bodies.
# stokes: Julia wrapper to call the Fortran stokessolver
function stokes(npts::Integer, nbods::Integer, xx::Vector{Float64}, yy::Vector{Float64})
	ntot = npts*nbods
	tau = zeros(Float64, ntot)
	ccall((:stokessolver_, "libstokes.so"),
		Void, (Ptr{Int}, Ptr{Int}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}), 
		&npts, &nbods, xx, yy, tau)
	return tau
end
# stokes!: Dispatch for vector of ThetaLenType; updates atau in each thlen.
function stokes!(thlenv::Vector{ThetaLenType})
	nbods = endof(thlenv)
	npts = endof(thlenv[1].theta)
	ntot = nbods*npts
	xv = zeros(Float64,ntot); yv = zeros(Float64,ntot)
	# Put all of the xy values in a single vector.
	for nn = 1:nbods
		# Compute the xy coordinates if they are not already loaded in thlen.
		getxy!(thlenv[nn])
		# Put the values into a single vector.
		n1 = npts*(nn-1)+1
		n2 = npts*nn
		xv[n1:n2], yv[n1:n2] = thlenv[nn].xx, thlenv[nn].yy
	end
	# Call the stokessolver.
	tau = stokes(npts,nbods,xv,yv)
	# Update the atau value in each of the thlen variables.
	for nn = 1:nbods
		n1 = npts*(nn-1)+1
		n2 = npts*nn
		thlenv[nn].atau = abs(tau[n1:n2])
	end
	return
end
##################################################

#################### Geometry functions ####################
# getalpha: Calculate the parameterization variable, alpha = s/L.
function getalpha(npts::Integer)
	# alpha = s/L is the parameterization variable.
	dalpha = 1.0/npts
	# Use an offset grid.
	alpha = collect(range(0.5*dalpha, dalpha, npts))
	return alpha
end
# circgeo: Creates a circle.
function circgeo(npts::Integer, rad::Float64, xsm::Float64=0.0, ysm::Float64=0.0)
	# Create a new ThetaLenType variable.
	thlen = new_thlen()
	# alpha = s/L is the parameterization variable.
	alpha = getalpha(npts)
	# theta is the tangent angle.
	thlen.theta = 0.5*pi + 2*pi*alpha
	# len is the total arclength.
	thlen.len = 2*pi*rad
	# Save xsm and ysm too.
	thlen.xsm = xsm; thlen.ysm = ysm
	return thlen
end
# polygongeo: Creates a polygon with number of sides, nsides.
function polygongeo(npts::Integer, nsides::Integer, 
		sigma::Float64 = 0.1, sdlen::Float64=0.5, xsm::Float64=0.0, ysm::Float64=0.0)
	# Create a new ThetaLenType variable.
	thlen = new_thlen()
	# alpha = s/L is the parameterization variable.
	alpha = getalpha(npts)
	# theta is the tangent angle.
	theta = zeros(npts)
	for nn=1:nsides
		a0 = (nn-1)/nsides
		a1 = nn/nsides
		intvl = ((a0 .<= alpha) & (alpha .<= a1))
		theta[intvl] = 0.5*pi + 2*pi*(nn-1)/nsides
	end
	# Smooth theta
	theta = gaussfilter(theta - 2*pi*alpha, sigma) + 2*pi*alpha
	# Save theta in the ThetaLenType object.
	thlen.theta = theta
	# len is the total arclength.
	thlen.len = nsides*sdlen
	# Save xsm and ysm too.
	thlen.xsm = xsm; thlen.ysm = ysm
	return thlen
end
##################################################
