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












function RK4(thld0::ThLenDenType, params::ParamType)
	
	epsilon = params.epsilon

	# Compute the derivatives k1.
	k1 = getderivs(thld0, params)
	# Compute the derivatives k2.
	thldtemp = feuler(thld0, 0.5*dt, k1, epsilon)
	k2 = getderivs(thldtemp, params)
	# Compute the derivatives k3.
	thldtemp = feuler(thld0, 0.5*dt, k2)
	k3 = getderivs(thldtemp, params)
	# Compute the derivatives k4.
	thldtemp = feuler(thld0, dt, k3)
	k4 = getderivs(thldtemp, params)
	# Compute the average derivatives and take the RK4 step.
	kavg = getkavg(k1,k2,k3,k4)
	thld1 = feuler(thld0, dt, kavg)
	return thld1
end
# feuler: Take a step of forward Euler with the derivatives specified by dvec.
function feuler(thld0::ThLenDenType, dt::Float64, dvec::Vector{DerivsType}, epsilon::Float64)
	npts, nbods = getnvals(thld0.thlenvec)
	thlv1 = new_thlenvec(nbods)
	for nn = 1:nbods
		# Extract the variables for body nn.
		thlen = thld0.thlenvec[nn]
		derivs = dvec[nn]
		th0, len0, xsm0, ysm0 = thlen.theta, thlen.len, thlen.xsm, thlen.ysm
		mterm, nterm = derivs.mterm, derivs.nterm
		xsmdot, ysmdot = derivs.xsmdot, derivs.ysmdot
		# Advance len with forward Euler first.
		len1 = len0 + dt*mterm
		# Advance theta with combination of integrating factor and forward Euler.
		sigma = 2*pi*sqrt(epsilon*dt*(elfun(len0)+elfun(len1)))
		th1 = gaussfilter(th0+dt*nterm, sigma)
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

#--------------- SMALL ROUTINES ---------------#
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
# getderivs: Get the vector of derivative terms for all the bodies.
function getderivs(thlenden::ThLenDenType, params::ParamType)
	getstress!(thlenden, params)
	npts, nbods = getnvals(thlenden.thlenvec)
	dvec = new_dvec(nbods)
	for nn = 1:nbods
		thlen = thlenden.thlenvec[nn]
		dvec[nn] = getderivs(thlen, params)
	end
	return dvec
end
#= getderivs: Calculate the derivatives of theta, len, xsm, ysm for a single body.
mterm=dL/dt; nterm is the nonlinear term in d/dt of theta. =#
function getderivs(thlen::ThetaLenType, params::ParamType)
	# Extract the variables.
	theta, len, atau = thlen.theta, thlen.len, thlen.atau
	epsilon = params.epsilon
	# Make sure atau has been computed.
	if atau==[]; throw("atau has not been computed"); return; end
	# Calculate the derivative terms.
	alpha = getalpha(endof(theta))
	dtheta = specdiff(theta - 2*pi*alpha) + 2*pi
	vnorm = atau + epsilon*cdfscale(len)*len^(-1) * (dtheta - 2*pi)
	vtang, mterm = tangvel(dtheta, vnorm)
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
	prod = vecmult(dtheta,vnorm)
	mprod = mean(prod)
	# Formula for mterm = dL/dt.
	mterm = -mprod
	# The derivative of the tangnential velocity wrt alpha.
	dvtang = prod - mprod
	# Spectrally integrate (mean-free) dvtan to get the tangential velocity.
	vtang = specint(dvtang)
	return vtang, mterm
end
# elfun: How to scale the smoothing with len.
function elfun(len::Float64)
	return 1./(len^2 * log(2*pi/len))
end











#--------------- MORE SMALL ROUTINES ---------------#
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
# trimthlenvec: Remove the curves with length too small or too big.
function trimthlenvec!(thlenvec1::Vector{ThetaLenType}, thlenvec0::Vector{ThetaLenType}, 
		minlen::Float64 = 1e-6, maxlen::Float64 = 2*pi)
	npts,nbods = getnvals(thlenvec1)
	lenvec = zeros(Float64,nbods)
	for nn=1:nbods
		lenvec[nn] = thlenvec1[nn].len
	end
	zind = find((lenvec.<minlen) | (lenvec.>maxlen))
	deleteat!(thlenvec0,zind)
	deleteat!(thlenvec1,zind)
end

# getns: Get the normal and tangent directions.
# Convention: CCW parameterization and inward pointing normal.
function getns(theta::Vector{Float64})
	# CCW tangent vector.
	sx = cos(theta)
	sy = sin(theta)
	# Inward pointing normal vector.
	nx = -sy
	ny = sx
	return sx,sy,nx,ny
end
# getns: Get the normal and tangent directions on the rotated grid.
function getnsrot(theta::Vector{Float64})
	# CCW tangent vector.
	sx = -sin(theta)
	sy = cos(theta)
	# Inward pointing normal vector.
	nx = -sy
	ny = sx
	return sx,sy,nx,ny
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

#--------------- De-aliasing multiplications ---------------#
#= This did not seem to improve stability any, so I voided it.
# mult_dealias: Multiply two vectors while upsampling to dealias.
function mult_dealias(uu::Vector{Float64}, vv::Vector{Float64})
	# Get the size of the vectors.
	npts = endof(uu)
	assert(endof(vv) == npts)
	# Get some integers.	
	nov2 = div(npts,2)
	nzeros = npts
	# Transform to spectral space.
	uhat = fft(uu)
	vhat = fft(vv)
	# Upsample by nzeros.
	uhat = [uhat[1:nov2]; zeros(nzeros); uhat[nov2+1:npts]]
	vhat = [vhat[1:nov2]; zeros(nzeros); vhat[nov2+1:npts]]
	# Return to physical space and multiply.
	uup, vup = ifft(uhat), ifft(vhat)
	imagtest(uup); imagtest(vup)
	uup, vup = real(uup), real(vup)
	product = uup .* vup
	# Transform to spectral space and downsample the product.
	prodhat = fft(product)
	deleteat!(prodhat, nov2+1:nov2+nzeros)
	# Transform back to physical space, test the imaginary part, and return.
	product = ifft(prodhat) * (npts+nzeros)/npts
	imagtest(product)
	return real(product)
end
=#
