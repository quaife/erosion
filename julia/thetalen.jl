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

#--------------- TIME-STEPPING ROUTINES ---------------#
# rungekutta4: Take a step forward with 4th order Runge-Kutta.
function rungekutta4(thld0::ThLenDenType, params::ParamType)
	# Extract parameters.
	# TO DO: specify dt adaptively
	dtmax = params.dt
	dt = dtmax
	epsilon = params.epsilon
	# Compute the derivatives k1.
	println("Stage 1 of Runge-Kutta")
	k1 = getderivs(thld0, params)
	# Check if the length is too small, if so, delete the body gracefully.
	checklen!(thld0, k1, dtmax)
	# Compute the derivatives k2.
	println("\nStage 2 of Runge-Kutta")
	thldtemp = feuler(thld0, 0.5*dt, k1, epsilon)
	k2 = getderivs(thldtemp, params)
	# Compute the derivatives k3.
	println("\nStage 3 of Runge-Kutta")
	thldtemp = feuler(thld0, 0.5*dt, k2, epsilon)
	k3 = getderivs(thldtemp, params)
	# Compute the derivatives k4.
	println("\nStage 4 of Runge-Kutta")
	thldtemp = feuler(thld0, dt, k3, epsilon)
	k4 = getderivs(thldtemp, params)
	# Compute the average derivatives and take the RK4 step.
	kavg = getkavg(k1,k2,k3,k4)
	thld1 = feuler(thld0, dt, kavg, epsilon)
	println("Completed Runge-Kutta step.")
	return thld1, dt
end
#= feuler: Take a step of forward Euler for all of the bodies with 
the time derivatives specified by dvec. =#
function feuler(thld0::ThLenDenType, dt::Float64, dvec::Vector{DerivsType}, epsilon::Float64)
	npts, nbods = getnvals(thld0.thlenvec)
	assert(endof(dvec) == nbods)
	thlv1 = new_thlenvec(nbods)
	for nn = 1:nbods
		# Extract the variables for body nn.
		thlen0 = thld0.thlenvec[nn]
		derivs = dvec[nn]
		th0, len0, xsm0, ysm0 = thlen0.theta, thlen0.len, thlen0.xsm, thlen0.ysm
		mterm, nterm = derivs.mterm, derivs.nterm
		xsmdot, ysmdot = derivs.xsmdot, derivs.ysmdot
		# Advance len with forward Euler first.
		len1 = len0 + dt*mterm
		# Make sure that the length is non-negative.
		assert(len1 >= 0.)
		# Advance theta with combination of integrating factor and forward Euler.
		sig1 = 2*pi*sqrt(epsilon*dt*(elfun(len0)+elfun(len1)))
		alpha = getalpha(npts)
		th1 = gaussfilter(th0-2*pi*alpha+dt*nterm, sig1) + 2*pi*alpha
		# Advance xsm and ysm with forward Euler.
		xsm1 = xsm0 + dt*xsmdot
		ysm1 = ysm0 + dt*ysmdot 
		# Save new values in thlenvec type.
		thlv1[nn].theta, thlv1[nn].len = th1, len1
		thlv1[nn].xsm, thlv1[nn].ysm = xsm1, ysm1
	end
	thld1 = new_thlenden(thlv1)
	return thld1
end
# getderivs: Get the derivative terms for all of the bodies.
function getderivs(thlenden::ThLenDenType, params::ParamType)
	getstress!(thlenden, params)
	nbods = endof(thlenden.thlenvec)
	dvec = new_dvec(nbods)
	for nn = 1:nbods
		thlen = thlenden.thlenvec[nn]
		dvec[nn] = getderivs(thlen, params.epsilon)
	end
	return dvec
end
#= getderivs: Get the derivative terms for a single body.
These are the derivatives of theta, len, xsm, and ysm.
mterm is dL/dt; nterm is the nonlinear term in dtheta/dt. =#
function getderivs(thlen::ThetaLenType, epsilon::Float64)
	# Extract the variables.
	theta, len, atau = thlen.theta, thlen.len, thlen.atau
	# Make sure atau has been computed.
	assert(endof(atau)>0)
	# Calculate the derivative terms.
	alpha = getalpha(endof(theta))
	dtheta = specdiff(theta - 2*pi*alpha) + 2*pi
	vnorm = atau + epsilon*len*elfun(len) * (dtheta - 2*pi)
	vtang, mterm = tangvel(dtheta, vnorm)
	# Check that mterm is negative.
	assert(mterm <= 0.)
	# Derivative of absolute-value of shear stress.
	datau = specdiff(atau)	
	nterm = (datau + vecmult(dtheta,vtang))/len
	# Get the derivatives of xsm and ysm.
	xsmdot = mean(-vecmult(vnorm,sin(theta)) + vecmult(vtang,cos(theta)))
	ysmdot = mean( vecmult(vnorm,cos(theta)) + vecmult(vtang,sin(theta)))
	derivs = DerivsType(mterm, nterm, xsmdot, ysmdot)
	return derivs
end
# tangvel: Compute the tangential velocity and mterm = dL/dt along the way.
function tangvel(dtheta::Vector{Float64}, vnorm::Vector{Float64})
	dthvn = vecmult(dtheta,vnorm)
	mdthvn = mean(dthvn)
	# Formula for mterm = dL/dt.
	mterm = -mdthvn
	# The derivative of the tangnential velocity wrt alpha.
	dvtang = dthvn - mdthvn
	# Spectrally integrate (mean-free) dvtan to get the tangential velocity.
	vtang = specint(dvtang)
	return vtang, mterm
end
# elfun: How to scale the smoothing with len.
function elfun(len::Float64)
	return 1./(len^2 * log(2*pi/len))
end
# function checklen!: Check if the length is too small, if so, delete the body gracefully.
function checklen!(thlenden::ThLenDenType, kvec::Vector{DerivsType}, dtmax::Float64)
	ndts = 1.	# The number of dt values to look into the future for len.
	nbods = endof(thlenden.thlenvec)
	assert(endof(kvec) == nbods)
	deletevec = Array(Int,0)
	for nn = 1:nbods
		len = thlenden.thlenvec[nn].len
		mterm = kvec[nn].mterm
		#= According to the scaling laws (neglecting the log term), len = -2*mterm*(tf-t) 
		So if len <= -2*mterm times some multiple of dtmax, the body will vanish soon. 
		If len is too small, delete the body in thlenden and its derivatives in kvec. =#
		if len <= -2*mterm*ndts*dtmax
			append!(deletevec,[nn])
		end
	end
	deleteat!(thlenden.thlenvec, deletevec)
	deleteat!(kvec, deletevec)
end
# getkavg: Average all of the k values for RK4
function getkavg(k1::Vector{DerivsType}, k2::Vector{DerivsType}, 
		k3::Vector{DerivsType}, k4::Vector{DerivsType})
	nbods = endof(k1)
	kavg = new_dvec(nbods)
	for nn=1:nbods
		mavg = (k1[nn].mterm + 2*k2[nn].mterm + 2*k3[nn].mterm + k4[nn].mterm)/6.
		navg = (k1[nn].nterm + 2*k2[nn].nterm + 2*k3[nn].nterm + k4[nn].nterm)/6.
		xdavg = (k1[nn].xsmdot + 2*k2[nn].xsmdot + 2*k3[nn].xsmdot + k4[nn].xsmdot)/6.
		ydavg = (k1[nn].ysmdot + 2*k2[nn].ysmdot + 2*k3[nn].ysmdot + k4[nn].ysmdot)/6.
		kavg[nn].mterm, kavg[nn].nterm = mavg, navg
		kavg[nn].xsmdot, kavg[nn].ysmdot = xdavg, ydavg 
	end
	return kavg
end

#--------------- SMALL ROUTINES ---------------#
# vecmult: Multiply two vectors with or without dealiasing.
function vecmult(uu::Vector{Float64},vv::Vector{Float64})
#	return mult_dealias(uu,vv)
	return uu.*vv
end
# getalpha: Calculate the parameterization variable, alpha = s/L, using an offset grid.
function getalpha(npts::Integer)
	dalpha = 1.0/npts
	return alpha = collect(range(0.5*dalpha, dalpha, npts))
end
# getxy!: Dispatch for input of type ThetaLenType. Only computes if they are not loaded.
function getxy!(thlen::ThetaLenType)
	if thlen.xx==[] || thlen.yy==[]
		thlen.xx, thlen.yy = getxy(thlen.theta, thlen.len, thlen.xsm, thlen.ysm)
	end
	return
end
#= getxy: Given theta and len, reconstruct the x and y coordinates of a body.
xsm and ysm are the boundary-averaged values. =#
function getxy(theta::Vector{Float64}, len::Float64, xsm::Float64, ysm::Float64)
	test_theta_means(theta)
	# The partial derivatives dx/dalpha and dy/dalpha.
	dx = len * (cos(theta) - mean(cos(theta)))
	dy = len * (sin(theta) - mean(sin(theta)))
	# Integrate to get the x,y coordinates; result will have mean zero.
	xx = specint(dx); yy = specint(dy)
	# Move to have the correct average values.
	xx += xsm; yy += ysm
	return xx,yy
end

#--------------- Tests for the theta vector ---------------#
#= test_theta_means: Test that cos(theta) and sin(theta) have zero mean.
These are conditions for theta to describe a closed curve 
in the equal arc length frame. =#
function test_theta_means(theta::Vector{Float64})
	npts = endof(theta)
	m1 = mean(cos(theta))
	m2 = mean(sin(theta))
	maxmean = maximum(abs([m1,m2]))
	thresh = 20./npts
	if maxmean > thresh
		warn("theta means")
		println("The max mean of sin, cos is: ", signif(maxmean,3), " > ", signif(thresh,3))
	end
	return
end
# test_theta_ends: Test that the difference between the first and last tangent angles is 2pi.
function test_theta_ends(theta::Vector{Float64}, thresh::Float64 = 0.2)
	# Use quadratic extrapolation to estimate theta at alpha=0 from both sides.
	th0left = 15/8*theta[1] - 5/4*theta[2] + 3/8*theta[3]
	th0right = 15/8*theta[end] - 5/4*theta[end-1] + 3/8*theta[end-2] - 2*pi
	# Compare the two extrpaolations.
	th0diff = abs(th0left - th0right)
	if th0diff > thresh
		throw("theta ends") 
		println("The difference between the ends is: ", signif(th0diff,3), " > ", signif(thresh,3))
	end
	return
end

