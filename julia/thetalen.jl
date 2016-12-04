#= thetalen.jl: Collection of codes for theta-L boundary evolution.

Variables
theta (vector): The tangent angle as it varies along the boundary.
len (scalar): The length of the curve, L.
atau (vector): The absolute value of the shear stress, computed by fluid solver.

Parameters
dt: The time-step.
epsilon: The constant for the smoothing, curvature-driven-flow component.
beta: The power of L hitting the curvature-term in the V_n; should be zero for Stokes flow. 

Convention for tangential and normal vectors
I use the same convention as Shelley 1994. That is, I assume the curve 
is parameterized in the counter-clockwise (CCW) direction, and I use the 
inward pointing normal. =#

#################### Multistep routines ####################
# advance_thetalen!: Dispatch for vectors of ThetaLenType to handle multiple bodies.
function advance_thetalen!(thlenvec1::Vector{ThetaLenType}, thlenvec0::Vector{ThetaLenType}, params::ParamType)
	# Call advance_thetalen for each element of the thlenvec.
	for ii = 1:endof(thlenvec0)
		advance_thetalen!(thlenvec1[ii],thlenvec0[ii],params)
	end
	# Remove curves with non-positive length.
	trimthlenvec!(thlenvec1, thlenvec0)
end
#= advance_thetalen: Advance theta and len from time-step n=1 to n=2.
This is a multi-step method and uses some values from n=0 too.
Note: thlen1 must already be loaded with the correct atua, which is used in getmn!() =#
function advance_thetalen!(thlen1::ThetaLenType, thlen0::ThetaLenType, params::ParamType)
	# Extract the needed variables.
	dt = params.dt; m0 = thlen0.mterm; len1 = thlen1.len
	# Calculate mterm and nterm at time n=1.
	m1 = getmn!(thlen1,params)
	# Update len with an explicit, multistep method; error if len negative.
	len2 = len1 + 0.5*dt*(3*m1-m0)
	# Create a new ThetaLenType variable and save the new len.
	thlen2 = new_thlen()
	thlen2.len = len2
	# Update theta with a multistep, integrating-factor method.
	advance_theta!(thlen2,thlen1,thlen0,params)
	# Update surface-mean coordinates with an explicit, multistep method.
	thlen2.xsm = thlen1.xsm + 0.5*dt*(3*thlen1.xsmdot - thlen0.xsmdot)
	thlen2.ysm = thlen1.ysm + 0.5*dt*(3*thlen1.ysmdot - thlen0.ysmdot)
	# Now thlen1 becomes the new thlen0, and thlen2 becomes the new thlen1
	copy_thlen!(thlen1,thlen0)
	copy_thlen!(thlen2,thlen1)
	return
end

# advance_theta: Advance theta in time with the integrating-factor method.
function advance_theta!(thlen2::ThetaLenType, thlen1::ThetaLenType, thlen0::ThetaLenType, params::ParamType)
	# Extract the needed variables.
	dt, epsilon, beta = params.dt, params.epsilon, params.beta
	theta1, len1, n1 = thlen1.theta, thlen1.len, thlen1.nterm
	len0, n0 = thlen0.len, thlen0.nterm
	len2 = thlen2.len
	alpha = getalpha(endof(theta1))
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
function getmn(theta::Vector{Float64}, len::Float64, atau::Vector{Float64}, params::ParamType)
	# Make sure that atau is not empty.
	if atau==[]; error("atau has not been computed"); return; end
	# Extract the needed variables.
	epsilon, beta = params.epsilon, params.beta
	alpha = getalpha(endof(theta))
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
	# Calculate the motion of the surface-mean points.
	xsmdot = mean(-vnorm.*sin(theta) + vtang.*cos(theta))
	ysmdot = mean( vnorm.*cos(theta) + vtang.*sin(theta))
	# Return everything.
	return mterm, nterm, xsmdot, ysmdot
end
#= getmn!: Dispatch for ThetaLenType input; saves mterm and nterm in thlen.
Also returms mterm to be used locally.
Note: thlen must already be loaded with the correct atau. =#
function getmn!(thlen::ThetaLenType, params::ParamType)
	thlen.mterm, thlen.nterm, thlen.xsmdot, thlen.ysmdot = 
		getmn(thlen.theta, thlen.len, thlen.atau, params)
	return thlen.mterm
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
# thetadot: Calculate the time derivative of theta; only used in the starter routine.
function thetadot(theta::Vector{Float64}, len::Float64, atau::Vector{Float64}, params::ParamType)
	# Extract the needed variables.
	dt, epsilon, beta = params.dt, params.epsilon, params.beta
	alpha = getalpha(endof(theta))
	# Calculate mterm and nterm at time 0.
	mterm, nterm, xsmdot, ysmdot = getmn(theta,len,atau,params)
	# Calculate the time derivative of theta.
	dth = specdiff(theta - 2*pi*alpha) + 2*pi
	d2th = specdiff(dth)
	thdot = epsilon*len^(beta-2)*d2th + nterm
	# Return thdot, mterm, nterm.
	return thdot, mterm, nterm, xsmdot, ysmdot
end
#= thetadot: Dispatch for ThetaLenType.
It also calculates mterm, nterm, xsmdot, ysmdot and saves them in thlen. =#
function thetadot!(thlen::ThetaLenType, params::ParamType)
	# Call thetadot and save mterm and nterm in thlen.
	thdot, thlen.mterm, thlen.nterm, thlen.xsmdot, thlen.ysmdot = 
		thetadot(thlen.theta, thlen.len, thlen.atau, params)
	return thdot
end
# trimthlenvec: Remove the curves with non-positive length.
function trimthlenvec!(thlenvec1::Vector{ThetaLenType}, thlenvec0::Vector{ThetaLenType})
	vsz = endof(thlenvec1)
	lenvec = zeros(Float64,vsz)
	for ii=1:vsz
		lenvec[ii] = thlenvec1[ii].len
	end
	zind = find(lenvec.<=0)
	deleteat!(thlenvec0,zind)
	deleteat!(thlenvec1,zind)
end