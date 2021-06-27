#-----------------------------------------------------------#
#= MAIN PURPOSE: Collection of routines for theta-L boundary evolution.

Variables
theta (vector): The tangent angle as it varies along the boundary of a given body.
len (scalar): The total arch length of the boundary of a given body.
xsm, ysm (scalars): The surface-mean of x and y; used as a reference point for each body.
atau (vector): The absolute value of the shear stress, computed by fluid solver.
matau (scalar): The surface mean of atau on each body.

Parameters
dt: The time-step.
epsilon: The constant for the smoothing, curvature-driven-flow component.
sigma: The smooth of the abolute-value of shear stress.

Parameterization Convention
The curve is assumed to be parametrized in the counter-clockwise (CCW) direction,
and the inward pointing normal vector is used. This is the same convectino as Shelley 1994.

Indexing Convention
bod = 1:nbods indexes the bodies. No other indices are used.

Exported Methods
rungekutta2() is the main time-stepping routine called in main.erosion().
 =#
#-----------------------------------------------------------#

# Create a module.
module TimeStepping
export rungekutta2

# Includes
using Statistics
include("spectral.jl")	# Used routines: specdiff, specint, expsmooth



# WILL BE MODULE?
include("callFortran.jl") # Used methods: ThetaLen Types, getstress!


#------------ ROUTINES TO COMPUTE TIME DERIVATIVES OF VARIABLES. ------------#
# Dot product of two vectors. If desired, I could de-aliase here.
vec_mult(uu::Vector{Float64}, vv::Vector{Float64}) = uu .* vv

# Compute the tangential velocity and mterm = dL/dt along the way.
function tangvel(dtheta::Vector{Float64}, vnorm::Vector{Float64})
	dthvn = vec_mult(dtheta, vnorm)
	mdthvn = mean(dthvn)
	# Formula for mterm = dL/dt.
	mterm = -2*pi*mdthvn
	# The derivative of the tangnential velocity wrt alpha.
	dvtang = dthvn .- mdthvn
	# Spectrally integrate (mean-free) dvtan to get the tangential velocity.
	vtang = specint(dvtang)
	return vtang, mterm
end

# Get alpha = 2 pi s/L, using an offset grid.
function getalpha(npts::Integer)
	dalpha = 2*pi/npts
	return alpha = collect(range(0.5*dalpha, step=dalpha, length=npts))
end

# Structure that includes the derivatives of theta, len, xsm, ysm.
mutable struct DerivsType
	mterm::Float64; nterm::Vector{Float64}; 
	xsmdot::Float64; ysmdot::Float64
end

#= Get the derivative terms for a single body.
These are the derivatives of theta, len, xsm, and ysm.
mterm is dL/dt; nterm is the nonlinear term in dtheta/dt. =#
function getderivs(thlen::ThetaLenType, params::ParamSet, atau::Vector{Float64})
	# Extract the variables.
	theta, len, matau = thlen.theta, thlen.len, thlen.matau
	@unpack epsilon, fixarea = params
	# Compute the effective atau to be used in evolution, depending on fixarea.
	ataueff = atau .- fixarea*matau
	# Calculate the derivative terms.
	alpha = getalpha(length(theta))
	dtheta = specdiff(theta - alpha) .+ 1.0
	vnorm = ataueff + epsilon*matau*(dtheta .- 1.0)
	vtang, mterm = tangvel(dtheta, vnorm)
	# Derivative of absolute-value of shear stress.
	datau = specdiff(atau)
	nterm = 2*pi/len * (datau + vec_mult(dtheta,vtang))
	# Get the derivatives of xsm and ysm.
	xsmdot = mean(-vec_mult(vnorm, sin.(theta)) + vec_mult(vtang, cos.(theta)))
	ysmdot = mean( vec_mult(vnorm, cos.(theta)) + vec_mult(vtang, sin.(theta)))
	# Save in object.
	derivs = DerivsType(mterm, nterm, xsmdot, ysmdot)
	return derivs
end

# Get the derivative terms for all of the bodies.
function getderivs(thlenden::ThLenDenType, params::ParamSet)
	atau = getstress!(thlenden, params)
	dvec = Array{DerivsType}(undef, 0)
	for bod = 1:length(thlenden.thlenvec)
		thlen = thlenden.thlenvec[bod]
		derivs = getderivs(thlen, params, atau[:,bod])
		push!(dvec, derivs)
	end
	return dvec
end
#-----------------------------------------------------------#


#--------------- TIME-STEPPING ROUTINES ---------------#
#= Find the bodies where the length is too small.
Neglecting the log term, L should vanish like sqrt(t).
The midpoint rule in RK2 gives the sqrt(2) factor. =#
function delete_indices(thld0::ThLenDenType, dvec::Vector{DerivsType}, dt::Float64)
	# ndts: The number of dt values to look into the future for len.
	ndts = 1.0
	mthresh = 10.
	thlv = thld0.thlenvec
	nbods = length(thlv)
	@assert length(dvec) == nbods
	deletevec = Array{Int}(undef,0)
	for bod = 1:nbods
		len = thlv[bod].len
		mterm = dvec[bod].mterm
		minlen = -sqrt(2)*mterm*ndts*dt
		if mterm > 0.; @warn("mterm is positive."); end;
		if (len <= minlen || mterm > mthresh)
			println("\n\n--------------------------------------------------")
			println("DELETING BODY ", bod)
			println("mterm = ", round(mterm,sigdigits=3), "; len = ", 
				round(len, sigdigits=3), "; minlen = ", round(minlen, sigdigits=3))
			println("--------------------------------------------------\n")
			append!(deletevec, [bod])
		end
	end
	return deletevec
end

#= Take a step of forward for all of the bodies.
dt1 is the step-size, dt2 is used in the Gaussian filter.  =#
function timestep!(thld0::ThLenDenType, thld_derivs::ThLenDenType, 
		dt1::Float64, dt2::Float64, params::ParamSet)
	# Compute the time derivatives.
	dvec = getderivs(thld_derivs, params)
	# Remove small bodies if needed.
	deletevec = delete_indices(thld0, dvec, dt1)
	deleteat!(thld0.thlenvec, deletevec)
	deleteat!(dvec, deletevec)
	# Only delete thld_derivs if it is different from thld0.
	if (thld0 != thld_derivs)
		deleteat!(thld_derivs.thlenvec, deletevec)
	end
	# Loop over all bodies and advance each forward in time.
	nbods = length(thld0.thlenvec)
	@assert (length(dvec) == length(thld_derivs.thlenvec) == nbods)
	thlv1 = Array{ThetaLenType}(undef,0)
	epsilon = params.epsilon
	alpha = getalpha(params.npts)
	# zetafun: How to scale the smoothing with len and matau = mean(abs(tau)).
	zetafun(len::Float64, matau::Float64) = 2*pi/len * matau
	for bod = 1:nbods
		# Extract the variables for the body of interest.
		thlen0  = thld0.thlenvec[bod]
		th0, len0, xsm0, ysm0 = thlen0.theta, thlen0.len, thlen0.xsm, thlen0.ysm
		# thlend	
		thlend = thld_derivs.thlenvec[bod]
		lend = thlend.len
		# matau = mean(atau)
		matau0, mataud = thlen0.matau, thlend.matau
		# derivs
		derivs = dvec[bod]
		mterm, nterm, xsmdot, ysmdot = derivs.mterm, derivs.nterm, derivs.xsmdot, derivs.ysmdot
		# Advance len first.
		len1 = len0 + dt1*mterm
		@assert (len1 > 0.)	
		# To advance theta, need zeta.
		zeta0 = zetafun(len0,matau0)
		zetad = zetafun(lend,mataud)
		# The first factor in the exponential smoothing.
		fac1 = epsilon*dt1*zetad
		# The second factor, for which there are two possible approaches.
		#fac2 = epsilon*dt2*0.5*(3*zetad-zeta0)	# Old Rule: uses trapezoid and midpoint combination.
		fac2 = epsilon*dt2*(2*zetad-zeta0)		# New Rule: uses RK2 for everything.
		# Advance to get the next theta.
		th1 = expsmooth(th0-alpha, fac1) + alpha
		th1 += dt1*expsmooth(nterm, fac2)
		# Advance xsm and ysm with forward Euler.
		xsm1 = xsm0 + dt1*xsmdot
		ysm1 = ysm0 + dt1*ysmdot 
		# Save the new thlen values in a vector.
		thlen1 = ThetaLenType(th1, len1, xsm1, ysm1, NaN)
		push!(thlv1,thlen1)
	end
	thld1 = new_thlenden(thlv1)
	return thld1
end

#= Take a step forward with 4th order Runge-Kutta.
This is the main time-stepping routine called by the erosion simulation. =#
function rungekutta2(thld0::ThLenDenType, params::ParamSet)
	# Could in principle set dt adaptively
	dt = params.dt
	# Stage 1
	println("Stage 1 of Runge-Kutta")
	thld05 = timestep!(thld0, thld0, 0.5*dt, 0.5*dt, params)
	# Stage 2
	println("\nStage 2 of Runge-Kutta")	
	thld1 = timestep!(thld0, thld05, dt, 0.5*dt, params)
	println("Completed Runge-Kutta step.")
	# Update the time value
	thld1.tt = thld0.tt + dt
	return thld1
end
#-----------------------------------------------------------#