#= thetalen.jl: Collection of codes for theta-L boundary evolution.
All routines work for only a single body unless stated otherwise.

Variables
atau (vector): The absolute value of the shear stress, computed by fluid solver.
theta (vector): The tangent angle as it varies along the boundary.
len (scalar): The length of the curve, L.
alpha(vector): Normalized arclength, s/L

Parameters
dt: The time-step.
epsilon: The constant for the smoothing, curvature-driven-flow component.
beta: The power of L hitting the curvature-term in the V_n; should be zero for Stokes flow. 

Convention for tangential and normal vectors
I use the same convention as Shelley 1994. That is, I assume the curve 
is parameterized in the counter-clockwise (CCW) direction, and I use the 
inward pointing normal. =#

#= advance_thetalen: Advance theta and len from time-step n=1 to n=2, 
using some values from n=0. This is a multistep method.
Note: thlen1 must already be loaded with the correct atua, which is used in getmn!() =#
function advance_thetalen!(thlen1::ThetaLenType, thlen0::ThetaLenType, params::ParamType)
	# Extract the needed variables.
	dt = params.dt
	m0 = thlen0.mterm
	len1 = thlen1.len
	# Calculate mterm and nterm at time n=1.
	getmn!(thlen1,params)
	m1 = thlen1.mterm
	# Update len with an explicit, multistep method; error if len negative.
	len2 = len1 + 0.5*dt*(3*m1-m0)
	if len2<0; error("The curve length is negative"); return; end
	# Create a new ThetaLenType variable and save the new len.
	thlen2 = new_thlen(thlen1)
	thlen2.len = len2
	# Update theta with a multistep, integrating-factor method.
	advance_theta!(thlen2,thlen1,thlen0,params)
	# Now thlen1 becomes the new thlen0, and thlen2 becomes the new thlen1
	thlen0 = deepcopy(thlen1)
	thlen1 = deepcopy(thlen2)
	return
end

# advance_theta: Advance theta in time with the integrating-factor method.
function advance_theta!(thlen2::ThetaLenType, thlen1::ThetaLenType, thlen0::ThetaLenType, params::ParamType)
	# Extract the needed variables.
	dt, epsilon, beta = params.dt, params.epsilon, params.beta
	theta1, len1, n1 = thlen1.theta, thlen1.len, thlen1.nterm
	len0, n0 = thlen0.len, thlen0.nterm
	len2 = thlen2.len
	alpha = thlen0.alpha
	# The power of L that is used.
	lpow = beta-2
	# The first value used in the Gaussian filter.
	sig1 = sqrt( len1^lpow + len2^lpow )
	sig1 *= 2*pi*sqrt(epsilon*dt)
	# The second value used in the Gaussian filter.
	sig2 = sqrt( len0^lpow + 2*len1^lpow + len2^lpow )
	sig2 *= 2*pi*sqrt(epsilon*dt)
	# Apply the appropriate Gaussian filters to advance theta in time.
	thlen2.theta = gaussfilter( theta1 - 2*pi*alpha, sig1) + 2*pi*alpha
	thlen2.theta += 0.5*dt*( 3*gaussfilter(n1,sig1) - gaussfilter(n0,sig2) )
	return
end

# getmn: Calculates mterm and nterm: mterm=dL/dt and nterm is the nonlinear term.
# Note: thlen must already be loaded with the correct atau.
function getmn!(thlen::ThetaLenType, params::ParamType)
	# Extract the needed variables.
	epsilon, beta = params.epsilon, params.beta
	alpha, theta, len, atau = thlen.alpha, thlen.theta, thlen.len, thlen.atau
	# The derivative of theta wrt alpha.
	dtheta = specdiff(theta - 2*pi*alpha) + 2*pi
	# The normal velocity.
	vnorm = atau + epsilon*len^(beta-1) * (dtheta - 2*pi)
	# Get the tangential velocity and dL/dt.
	vtang, mterm = tangvel(dtheta, vnorm)
	# The derivative of the absolute value of shear stress.
	datau = specdiff(atau)
	# The nonlinear term in the theta-evolution equation.
	nterm = (datau + dtheta.*vtang)/len
	# Save the results in the thlen variable.
	thlen.mterm = mterm
	thlen.nterm = nterm
	return
end

# tangvel: Compute the tangential velocity and mterm = dL/dt along the way.
function tangvel(dtheta::Vector{Float64}, vnorm::Vector{Float64})
	# The product of dtheta and vnorm, along with its mean.
	prod = dtheta.*vnorm
	mprod = mean(prod)
	# Formula for mterm = dL/dt.
	mterm = -mprod
	# The derivative of the tangnential velocity wrt alpha.
	dvtang = prod - mprod
	# Spectrally integrate (mean-free) dvtan to get the tangential velocity.
	vtang = specint(dvtang)
	return vtang, mterm
end
