# OBJECTIVE: General use objects and methods related to theta-len variables.
# Convention: bod = 1:nbods indexes the bodies.

module ThetaLen
export ParamSet, ThetaLenType, ThLenDenType, new_thlenden, getxy, getnxy

using ..SpectralMethods: specint

#-------------------------------------------------#

using Statistics: mean
using Parameters: @with_kw

#--------------- ParamSet data type ---------------#
# Collect all of the parameters with sane defaults.
@with_kw struct ParamSet
	# The data file for the initial geometry
	infolder::AbstractString = "input_geos/test/"
	label::AbstractString = "02-1"
	# Varied computational parameters
	npts::Int = 128			# The number of points per body, default 128.
	ibary::Int = 1			# Use barycentric (1) or not (0).
	ifmm::Int = 1			# Use the FFM (1) or not (0).
	ibc::Int = 1			# Use slip BCs (1) or no-slip (0)
	# Varied physical parameters
	epsfac::Float64 = 15	# Smoothing parameter for curvature driven flow.
	sigfac::Float64 = 10	# Smoothing parameter for the stress.
	dt::Float64 = 1e-3		# The time step.
	outstride::Int = 4		# The number of steps per output, default 4.
	# Fixed physical parameters
	fixpdrop::Bool = 1		# Fix the pressure drop (1) or not (0)
	fixarea::Bool = 0		# Keep the area fixed (1) or not (0)
	tfin::Float64 = 1.0		# The final time.
	# Fixed computational parameters
	maxl::Int = 8000		# The maximum number of GMRES iterations.
	nouter::Int = 1024		# The number of points on the outer boundary, default 1024.
	# Derived parameters
	epsilon::Float64 = epsfac/npts
	sigma::Float64 = sigfac/npts
end

#--------------- Theta-Len data types ---------------#
#= Structure that includes the geometry variables (theta, len, xsm, ysm),
and matau: the mean of the absolute stress on each body. =#
mutable struct ThetaLenType
	theta::Vector{Float64}; len::Float64; xsm::Float64; ysm::Float64; matau::Float64
end
# Structure that includes the vector of thlens for all bodies and the density function.
mutable struct ThLenDenType
	thlenvec::Vector{ThetaLenType}; tt::Float64;
	density::Vector{Float64}; denrot::Vector{Float64}; 
end
# Create a new ThLenDenType variable.
new_thlenden(thlenvec::Vector{ThetaLenType}) = ThLenDenType(thlenvec, 0., [],[])

#--------------- Routines for x-y coordinates ---------------#
# Given theta and len, reconstruct the x-y coordinates of a single body.
# xsm, ysm are the surface-mean values used as a reference point.
function getxy(thlen::ThetaLenType)
	theta, len, xsm, ysm = thlen.theta, thlen.len, thlen.xsm, thlen.ysm
	test_theta_means(theta)
	@assert len > 0.
	dx = len/(2*pi) * (cos.(theta) .- mean(cos.(theta)))
	dy = len/(2*pi) * (sin.(theta) .- mean(sin.(theta)))
	xx = specint(dx) .+ xsm
	yy = specint(dy) .+ ysm
	return xx, yy
end

# Get the x-y coordinates of all the bodies and also nbods.
function getnxy(thlenden::ThLenDenType)
	nbods = length(thlenden.thlenvec)
	xv,yv = [Array{Float64}(undef,0) for ii=1:2]
	for bod = 1:nbods
		xx, yy = getxy(thlenden.thlenvec[bod])
		append!(xv,xx); append!(yv,yy)
	end
	return nbods, xv, yv
end

#--------------- Tests for the theta vector ---------------#
# These test appear to be UNUSED, except possibly in remake_data.
#= Test that cos(theta) and sin(theta) have zero mean.
These are conditions for theta to describe a closed curve in the equal arc length frame. =#
function test_theta_means(theta::Vector{Float64})
	maxmean = maximum(abs, [mean(cos.(theta)), mean(sin.(theta))])
	thresh = 20.0 / length(theta)
	if maxmean > thresh
		@warn("theta means")
		println("The maximum mean of sin(theta) & cos(theta) is: ", 
			round(maxmean, sigdigits=3), " > ", round(thresh, sigdigits=3))
	end
end
#= Test that the difference between the first and last tangent angles is 2*pi. 
Current, this routine is UNUSED. =#
function test_theta_ends(theta::Vector{Float64}, thresh::Float64 = 0.2)
	# Use quadratic extrapolation to estimate theta at alpha=0 from both sides.
	th0left = 15/8*theta[1] - 5/4*theta[2] + 3/8*theta[3]
	th0right = 15/8*theta[end] - 5/4*theta[end-1] + 3/8*theta[end-2] - 2*pi
	# Compare the two extrpaolations.
	th0diff = abs(th0left - th0right)
	if th0diff > thresh
		@warn("theta ends") 
		println("The difference between the ends is: ", 
			round(th0diff, sigdigits=3), " > ", round(thresh, sigdigits=3))
	end
end
#-----------------------------------------------------------------------#

end