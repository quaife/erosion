#= thetalen.jl: Collection of codes for theta-L boundary evolution.

Variables
theta (vector): The tangent angle as it varies along the boundary.
len (scalar): The length of the curve, L.
atau (vector): The absolute value of the shear stress, computed by fluid solver.

Parameters
dt: The time-step.
epsilon: The constant for the smoothing, curvature-driven-flow component.
sigma: The smooth of the abolute-value of shear stress.

I use the same convention for tangential and normal vectors as Shelley 1994. 
That is, I assume the curve is parameterized in the counter-clockwise (CCW) direction, 
and I use the inward pointing normal vector. =#

#= cdfscale: The function to scale the curvature-driven flow appropriately with the shear stress. 
2D Stokes: -1/log(L); 3D Stokes: 1; high Reynolds: sqrt(L) 
Note: for 2D Stokes, I have to use log-of-tanh to get the same behvaior near L=0, 
but avoid problems at len=1, log(1) = 0. =#
function cdfscale(len::Float64)
	return -1./log(0.5*tanh(2*len))
end

#################### Multistep routines ####################
# advance_thetalen!: Dispatch for ThLenDenType.
function advance_thetalen!(thlenden1::ThLenDenType, thlenden0::ThLenDenType, params::ParamType)
	advance_thetalen!(thlenden1.thlenvec, thlenden0.thlenvec)
end
# advance_thetalen!: Dispatch for vectors of ThetaLenType to handle multiple bodies.
function advance_thetalen!(thlenvec1::Vector{ThetaLenType}, thlenvec0::Vector{ThetaLenType}, params::ParamType)
	# Call advance_thetalen for each element of the thlenvec.
	for nn = 1:endof(thlenvec0)
		advance_thetalen!(thlenvec1[nn],thlenvec0[nn],params)
	end
	# Remove curves with non-positive length.
	trimthlenvec!(thlenvec1, thlenvec0)
end
#= advance_thetalen: Advance theta and len from time-step n=1 to n=2.
This is a multi-step method and uses some values from n=0 too.
Note: thlen1 must already be loaded with the correct atua, which is used in getmn!() =#
function advance_thetalen!(thlen1::ThetaLenType, thlen0::ThetaLenType, params::ParamType)
	# Extract the needed variables.
	dt, m0, len1 = params.dt, thlen0.mterm, thlen1.len
	# Calculate mterm and nterm at time n=1.
	m1 = getmn!(thlen1,params.epsilon)
	# Update len with an explicit, multistep method.
	len2 = len1 + 0.5*dt*(3*m1-m0)
	# Create a new ThetaLenType variable and save the new len.
	thlen2 = new_thlen()
	thlen2.len = len2
	# Update theta with a multistep, integrating-factor method.
	advance_theta!(thlen2,thlen1,thlen0,params)
	# Update surface-mean coordinates with an explicit, multistep method.
	thlen2.xsm = thlen1.xsm + 0.5*dt*(3*thlen1.xsmdot - thlen0.xsmdot)
	thlen2.ysm = thlen1.ysm + 0.5*dt*(3*thlen1.ysmdot - thlen0.ysmdot)
	# Now thlen1 becomes the new thlen0, and thlen2 becomes the new thlen1.
	thlen0 = deepcopy(thlen1)
	thlen1 = deepcopy(thlen2)
	return
end

# advance_theta: Advance theta in time with the integrating-factor method.
function advance_theta!(thlen2::ThetaLenType, thlen1::ThetaLenType, thlen0::ThetaLenType, params::ParamType)
	# Extract the needed variables.
	dt, epsilon = params.dt, params.epsilon
	theta1, len1, n1 = thlen1.theta, thlen1.len, thlen1.nterm
	len0, n0 = thlen0.len, thlen0.nterm
	len2 = thlen2.len
	alpha = getalpha(endof(theta1))
	# The function that enters the sigmas of the Guassian filter.
	lenfun(len::Float64) = len^(-2)*cdfscale(len)
	# The first value used in the Gaussian filter.
	sig1 = sqrt( lenfun(len1) + lenfun(len2) )
	sig1 *= 2*pi*sqrt(epsilon*dt)
	# The second value used in the Gaussian filter.
	sig2 = sqrt( lenfun(len0) + 2*lenfun(len1) + lenfun(len2) )
	sig2 *= 2*pi*sqrt(epsilon*dt)
	# Apply the appropriate Gaussian filters to advance theta in time.
	thlen2.theta = gaussfilter( theta1 - 2*pi*alpha, sig1) + 2*pi*alpha
	thlen2.theta += 0.5*dt*( 3*gaussfilter(n1,sig1) - gaussfilter(n0,sig2) )
	return
end

#= getmn!: Dispatch for ThetaLenType input; saves mterm and nterm in thlen.
Also returms mterm to be used locally.
Note: thlen must already be loaded with the correct atau. =#
function getmn!(thlen::ThetaLenType, epsilon::Float64)
	thlen.mterm, thlen.nterm, thlen.xsmdot, thlen.ysmdot = 
		getmn(thlen.theta, thlen.len, thlen.atau, epsilon)
	return thlen.mterm
end
# getmn: Calculates mterm and nterm: mterm=dL/dt and nterm is the nonlinear term.
function getmn(theta::Vector{Float64}, len::Float64, atau::Vector{Float64}, epsilon::Float64)
	if atau==[]; throw("atau has not been computed"); return; end
	alpha = getalpha(endof(theta))
	dtheta = specdiff(theta - 2*pi*alpha) + 2*pi
	vnorm = atau + epsilon*cdfscale(len)*len^(-1) * (dtheta - 2*pi)
	vtang, mterm = tangvel(dtheta, vnorm)
	datau = specdiff(atau)	# Derivative of absolute-value of shear stress.
	nterm = (datau + dtheta.*vtang)/len
	xsmdot = mean(-vnorm.*sin(theta) + vtang.*cos(theta))
	ysmdot = mean( vnorm.*cos(theta) + vtang.*sin(theta))
	return mterm, nterm, xsmdot, ysmdot
end
# tangvel: Compute the tangential velocity and mterm = dL/dt along the way.
function tangvel(dtheta::Vector{Float64}, vnorm::Vector{Float64})
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
# trimthlenvec: Dispatch for ThLenDenType.
function trimthlenvec!(thlenden1::ThLenDenType, thlenden0::ThLenDenType)
	trimthlenvec!(thlenden1.thlenvec, thlenden0.thlenvec)
end
# trimthlenvec: Remove the curves with non-positive length.
function trimthlenvec!(thlenvec1::Vector{ThetaLenType}, thlenvec0::Vector{ThetaLenType})
	nbods = endof(thlenvec1)
	lenvec = zeros(Float64,nbods)
	for nn=1:nbods
		lenvec[nn] = thlenvec1[nn].len
	end
	zind = find(lenvec.<=0)
	deleteat!(thlenvec0,zind)
	deleteat!(thlenvec1,zind)
end
