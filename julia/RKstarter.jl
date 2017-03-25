# RKstarter.jl: Runge-Kutta Starter routines.

#= RKstarter!: Explicit second-order Runge-Kutta to start the time stepping.
It also calculates mterm, nterm, xsmdot, ysmdot and saves them in thlenden0. =#
function RKstarter!(thlenden0::ThLenDenType, params::ParamType)
	dt = params.dt
	epsilon = params.epsilon

	println("stop B, ", endof(thlenden0.thlenvec))

	# Compute the stress at t=0 and take the first step of RK2.
	getstress!(thlenden0, params)

	println("stop C, ", endof(thlenden0.thlenvec))

	thlenden05 = festep(0.5*dt, thlenden0, thlenden0, epsilon)
	# Compute the stress at t=0.5*dt and take the second step of RK2.
	getstress!(thlenden05, params)
	thlenden1 = festep(dt, thlenden0, thlenden05, epsilon)
	# Remove curves with non-positive length and return.
	trimthlenvec!(thlenden1, thlenden0)
	return thlenden1
end

# festep: Dispatch for ThLenDenType.
# Take an FE-step for each component of thlenvec.
function festep(dt::Float64, thlenden0::ThLenDenType, thlendendots::ThLenDenType, epsilon::Float64)
	nbods = endof(thlenden0.thlenvec)
	thlenden1 = new_thlenden()
	for nn = 1:nbods
		thlen0 = thlenden0.thlenvec[nn]
		thlendots = thlendendots.thlenvec[nn]
		thlenden1.thlenvec[nn] = festep(dt, thlen0, thlendots, epsilon)
	end
	return thlenden1
end
# festep: Dispatch for ThetaLenType.
# Allow the starting point and the point where the derivatives are taken to be different.
function festep(dt::Float64, thlen0::ThetaLenType, thlendots::ThetaLenType, epsilon::Float64)
	thdot = thetadot!(thlendots, epsilon)
	thlen1 = new_thlen()
	thlen1.theta, thlen1.len, thlen1.xsm, thlen1.ysm = festep(dt, thdot, 
		thlen0.theta, thlen0.len, thlen0.xsm, thlen0.ysm, 
		thlendots.mterm, thlendots.xsmdot, thlendots.ysmdot)
	return thlen1
end
#= festep: Take a single forward Euler step of theta, len, xsm, and ysm.
Note: Do not use the dt inside params because I might want to input 
something else, like 0.5*dt for the Runge-Kutta starter. =#
function festep(dt::Float64, thdot::Vector{Float64}, 
		theta0::Vector{Float64}, len0::Float64, xsm0::Float64, ysm0::Float64,
		ldot::Float64, xsmdot::Float64, ysmdot::Float64)
	theta1 = theta0 + dt*thdot
	len1 = len0 + dt*ldot
	xsm1 = xsm0 + dt*xsmdot
	ysm1 = ysm0 + dt*ysmdot
	return theta1, len1, xsm1, ysm1
end
#= thetadot: Dispatch for ThetaLenType.
It also calculates mterm, nterm, xsmdot, ysmdot and saves them in thlen. =#
function thetadot!(thlen::ThetaLenType, epsilon::Float64)
	thdot, thlen.mterm, thlen.nterm, thlen.xsmdot, thlen.ysmdot = 
		thetadot(thlen.theta, thlen.len, thlen.atau, epsilon)
	return thdot
end
# thetadot: Calculate the time derivative of theta; only used in the starter routine.
function thetadot(theta::Vector{Float64}, len::Float64, atau::Vector{Float64}, epsilon::Float64)
	# Calculate mterm and nterm at time 0.
	mterm, nterm, xsmdot, ysmdot = getmn(theta,len,atau,epsilon)
	# Calculate the time derivative of theta.
	alpha = getalpha(endof(theta))
	dth = specdiff(theta - 2*pi*alpha) + 2*pi
	d2th = specdiff(dth)
	thdot = epsilon*cdfscale(len)*len^(-2)*d2th + nterm
	return thdot, mterm, nterm, xsmdot, ysmdot
end
