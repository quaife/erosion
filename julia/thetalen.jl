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


########## Multi-step routines ##########
#= advance_thetalen: Advance theta and len from time-step n=1 to n=2, 
using some values from n=0. This is a multistep method. =#
function advance_thetalen(atau1::Vector{Float64}, theta1::Vector{Float64}, 
		len0::Float64, len1::Float64, M0::Float64, N0::Vector{Float64}, 
		dt::Float64, epsilon::Float64, beta::Real)
	# Get the terms M and N at time n=1.
	M1,N1 = getMN(atau1,theta1,len1,epsilon,beta)
	# Update len with an explicit, multistep method.
	len2 = len1 + 0.5*dt*(3*M1-M0)
	if len2<0; error("The curve length is negative"); return; end
	# Update theta with a multistep, integrating-factor method.
	theta2 = advance_theta(theta1,len0,len1,len2,N0,N1,dt,epsilon,beta)
	return theta2, len1, len2, M1, N1
end

# advance_theta: Advance theta in time with the integrating-factor method.
function advance_theta(theta1::Vector{Float64}, len0::Float64, len1::Float64, len2::Float64, 
		N0::Vector{Float64}, N1::Vector{Float64}, dt::Float64, epsilon::Float64, beta::Real)
	# Calculate alpha.
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
	theta2 = gaussfilter( theta1 - 2*pi*alpha, sig1) + 2*pi*alpha
	theta2 += 0.5*dt*( 3*gaussfilter(N1,sig1) - gaussfilter(N0,sig2) )
	return theta2
end

# getMN: Calculates the terms M=dL/dt and N.
function getMN(atau::Vector{Float64}, theta::Vector{Float64}, len::Float64, 
		epsilon::Float64, beta::Real)
	# The derivative of theta wrt alpha.
	alpha = getalpha(endof(theta))
	dtheta = specdiff(theta - 2*pi*alpha) + 2*pi
	# Calculate the normal velocity.
	vnorm = atau + epsilon*len^(beta-1) * (dtheta - 2*pi)
	# Get the tangential velocity and dL/dt.
	vtang, dldt = tangvel(dtheta, vnorm)
	# The derivative of the absolute value of shear stress.
	datau = specdiff(atau)
	# Calculate the nonlinear term in the theta-evolution equation.
	nterm = (datau + dtheta.*vtang)/len
	return dldt, nterm
end

# tangvel: Compute the tangential velocity and dL/dt along the way.
function tangvel(dtheta::Vector{Float64}, vnorm::Vector{Float64})
	# The product of dtheta and vnorm, along with its mean.
	prod = dtheta.*vnorm
	mprod = mean(prod)
	# Formula for dL/dt.
	dldt = -mprod
	# The derivative of the tangnential velocity wrt alpha.
	dvtang = prod - mprod
	# Spectrally integrate (mean-free) dvtan to get the tangential velocity.
	vtang = specint(dvtang)
	return vtang, dldt
end


########## Starter routines ##########
#= thetadot: Calculate the time derivative of theta;
Also return MM and NN while we're at it. =#
function thetadot(atau::Vector{Float64}, theta::Vector{Float64}, len::Float64, 
		epsilon::Float64, beta::Real)
	# Get the M and N terms
	MM,NN = getMN(atau,theta,len,epsilon,beta)
	# Calculate the time derivative of theta.
	alpha = getalpha(endof(theta))
	dth = specdiff(theta - 2*pi*alpha) + 2*pi
	d2th = specdiff(dth)
	thdot = epsilon*len^(beta-2)*d2th + NN
	return thdot, MM, NN
end


########## Get x and y coordinates ##########
#= getxy: Given theta and len, reconstruct the x and y coordinates of a body.
xc and yc are the coordinates of the center of mass.
While we're at it, also calculate the normal direcations. =#
function getxy(theta::Vector{Float64}, len::Float64, 
		xc::Float64=0.0, yc::Float64=0.0)
	# The increments of dx and dy
	dx = len * (cos(theta) - mean(cos(theta)))
	dy = len * (sin(theta) - mean(sin(theta)))
	# Integrate to get the x,y coordinates.
	xx = specint(dx)
	yy = specint(dy)

	## TO DO: Move the center of mass or average values.

	# The normal vector, direction???
	nx = -sin(theta)
	ny = cos(theta)
	return xx,yy
end
