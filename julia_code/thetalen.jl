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
function rungekutta2(thld0::ThLenDenType, params::ParamType)
	# Extract parameters.
	dt = params.dt
	epsilon = params.epsilon
	# Stage 1
	println("Stage 1 of Runge-Kutta")
	thld05 = timestep!(thld0, thld0, 0.5*dt, 0.5*dt, params)
	# Stage 2
	println("\nStage 2 of Runge-Kutta")	
	thld1 = timestep!(thld0, thld05, dt, 0.5*dt, params)
	println("Completed Runge-Kutta step.")
	return thld1, dt
end
#= timestep: Take a step of forward for all of the bodies.
dt1 is the step-size, dt2 is used in the Gaussian filter.  =#
function timestep!(thld0::ThLenDenType, thld_derivs::ThLenDenType, 
		dt1::Float64, dt2::Float64, params::ParamType)
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
	npts, nbods = getnvals(thld0.thlenvec)
	@assert (length(dvec) == length(thld_derivs.thlenvec) == nbods)
	thlv1 = new_thlenvec(nbods)
	epsilon = params.epsilon
	alpha = getalpha(npts)
	for nn = 1:nbods
		# Extract the variables for body nn.
		# thlen0
		thlen0  = thld0.thlenvec[nn]
		th0, len0, xsm0, ysm0 = thlen0.theta, thlen0.len, thlen0.xsm, thlen0.ysm
		# derivs
		derivs = dvec[nn]
		mterm, nterm, xsmdot, ysmdot = derivs.mterm, derivs.nterm, derivs.xsmdot, derivs.ysmdot
		# thlend	
		thlend = thld_derivs.thlenvec[nn]
		lend = thlend.len
		# matau = mean(atau)
		matau0 = mean(thlen0.atau)
		mataud = mean(thlend.atau)	
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
		th1 = expsmooth(th0-alpha,fac1) + alpha
		th1 += dt1*expsmooth(nterm,fac2)
		# Advance xsm and ysm with forward Euler.
		xsm1 = xsm0 + dt1*xsmdot
		ysm1 = ysm0 + dt1*ysmdot 
		# Save new values in thlenvec type.
		thlv1[nn].theta, thlv1[nn].len = th1, len1
		thlv1[nn].xsm, thlv1[nn].ysm = xsm1, ysm1
	end
	thld1 = new_thlenden(thlv1)
	return thld1
end
# zetafun: How to scale the smoothing with len and matau = mean(abs(tau)).
function zetafun(len::Float64, matau::Float64)
	return 2*pi/len * matau
end

#--------------- ROUTINES TO SUPPORT TIMESTEPPING ---------------#
# getderivs: Get the derivative terms for all of the bodies.
function getderivs(thlenden::ThLenDenType, params::ParamType)
	getstress!(thlenden, params)
	nbods = length(thlenden.thlenvec)
	dvec = new_dvec(nbods)
	for nn = 1:nbods
		thlen = thlenden.thlenvec[nn]
		dvec[nn] = getderivs(thlen, params.epsilon, params.fixarea)
	end
	return dvec
end
#= getderivs: Get the derivative terms for a single body.
These are the derivatives of theta, len, xsm, and ysm.
mterm is dL/dt; nterm is the nonlinear term in dtheta/dt. =#
function getderivs(thlen::ThetaLenType, epsilon::Float64, fixarea::Bool)
	# Extract the variables.
	theta, len, atau = thlen.theta, thlen.len, thlen.atau
	# Make sure atau has been computed.
	@assert length(atau)>0
	# Compute the effective atau to be used in evolution, depending on fixarea.
	matau = mean(atau)
	ataueff = atau .- fixarea*matau
	# Calculate the derivative terms.
	alpha = getalpha(length(theta))
	dtheta = specdiff(theta - alpha) .+ 1.0
	vnorm = ataueff + epsilon*matau*(dtheta .- 1.0)
	vtang, mterm = tangvel(dtheta, vnorm)
	# Derivative of absolute-value of shear stress.
	datau = specdiff(atau)
	nterm = 2*pi/len * (datau + vecmult(dtheta,vtang))
	# Get the derivatives of xsm and ysm.
	xsmdot = mean(-vecmult(vnorm,sin.(theta)) + vecmult(vtang,cos.(theta)))
	ysmdot = mean( vecmult(vnorm,cos.(theta)) + vecmult(vtang,sin.(theta)))
	# Save in object.
	derivs = DerivsType(mterm, nterm, xsmdot, ysmdot)
	return derivs
end
# tangvel: Compute the tangential velocity and mterm = dL/dt along the way.
function tangvel(dtheta::Vector{Float64}, vnorm::Vector{Float64})
	dthvn = vecmult(dtheta,vnorm)
	mdthvn = mean(dthvn)
	# Formula for mterm = dL/dt.
	mterm = -2*pi*mdthvn
	# The derivative of the tangnential velocity wrt alpha.
	dvtang = dthvn .- mdthvn
	# Spectrally integrate (mean-free) dvtan to get the tangential velocity.
	vtang = specint(dvtang)
	return vtang, mterm
end
#= delete_indices: Find the bodies where the length is too small.
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
	for nn = 1:nbods
		len = thlv[nn].len
		mterm = dvec[nn].mterm
		minlen = -sqrt(2)*mterm*ndts*dt
		if mterm > 0.; @warn("mterm is positive."); end;
		if (len <= minlen || mterm > mthresh)
			println("\n\n--------------------------------------------------")
			println("DELETING BODY ", nn)
			println("mterm = ", round(mterm,sigdigits=3), "; len = ", 
				round(len,sigdigits=3), "; minlen = ", round(minlen,sigdigits=3))
			println("--------------------------------------------------\n")
			append!(deletevec,[nn])
		end
	end
	return deletevec
end

#--------------- SMALL ROUTINES ---------------#
# vecmult: Multiply two vectors with or without dealiasing.
function vecmult(uu::Vector{Float64},vv::Vector{Float64})
#	return mult_dealias(uu,vv)
	return uu.*vv
end
# getalpha: Get alpha = 2 pi s/L, using an offset grid.
function getalpha(npts::Integer)
	dalpha = 2*pi/npts
	return alpha = collect(range(0.5*dalpha, step=dalpha, length=npts))
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
	@assert len > 0.
	# The partial derivatives dx/dalpha and dy/dalpha.
	dx = len/(2*pi) * (cos.(theta) .- mean(cos.(theta)))
	dy = len/(2*pi) * (sin.(theta) .- mean(sin.(theta)))
	# Integrate to get the x,y coordinates; result will have mean zero.
	xx = specint(dx); yy = specint(dy)
	# Move to have the correct average values.
	xx .+= xsm; yy .+= ysm
	return xx,yy
end

#--------------- Tests for the theta vector ---------------#
#= test_theta_means: Test that cos(theta) and sin(theta) have zero mean.
These are conditions for theta to describe a closed curve 
in the equal arc length frame. =#
function test_theta_means(theta::Vector{Float64})
	npts = length(theta)
	m1 = mean(cos.(theta))
	m2 = mean(sin.(theta))
	maxmean = maximum(abs.([m1,m2]))
	thresh = 20/npts
	if maxmean > thresh
		@warn("theta means")
		println("The max mean of sin, cos is: ", 
			round(maxmean,sigdigits=3), " > ", round(thresh,sigdigits=3))
	end
	return
end
#= test_theta_ends: Test that the difference between the 
first and last tangent angles is 2pi. =#
function test_theta_ends(theta::Vector{Float64}, thresh::Float64 = 0.2)
	# Use quadratic extrapolation to estimate theta at alpha=0 from both sides.
	th0left = 15/8*theta[1] - 5/4*theta[2] + 3/8*theta[3]
	th0right = 15/8*theta[end] - 5/4*theta[end-1] + 3/8*theta[end-2] - 2*pi
	# Compare the two extrpaolations.
	th0diff = abs(th0left - th0right)
	if th0diff > thresh
		@warn("theta ends") 
		println("The difference between the ends is: ", 
			round(th0diff,sigdigits=3), " > ", round(thresh,sigdigits=3))
	end
	return
end