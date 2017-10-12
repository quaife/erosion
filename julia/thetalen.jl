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
2D Stokes: cdfscale = 1/log(W/L), where W > L is lengthscale of 3rd dimension.
3D Stokes, cdfscale = 1; high Reynolds, cdfscale = sqrt(L) =#
function cdfscale(len::Float64)
	return 1./log(2*pi/len)
end

#--------------- MULTISTEP METHOD FOR THETA AND LEN ---------------#
# advance_thetalen!: Dispatch for ThLenDenType.
function advance_thetalen!(thlenden1::ThLenDenType, thlenden0::ThLenDenType, params::ParamType)
	advance_thetalen!(thlenden1.thlenvec, thlenden0.thlenvec, params)
	thlenden0.density = thlenden1.density
	thlenden1.density = evec()
	thlenden1.denrot = evec()
	return
end
# advance_thetalen!: Dispatch for vectors of ThetaLenType to handle multiple bodies.
function advance_thetalen!(thlenvec1::Vector{ThetaLenType}, thlenvec0::Vector{ThetaLenType}, params::ParamType)
	# Call advance_thetalen for each element of the thlenvec.
	for nn = 1:endof(thlenvec0)
		advance_thetalen!(thlenvec1[nn],thlenvec0[nn],params)
	end
	# Remove curves with non-positive length.
	trimthlenvec!(thlenvec1, thlenvec0)
	return
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
	# Now thlen1 gets copied to thlen0, and thlen2 gets copied to thlen1.
	# Note: Although tempting, I cannot use deepcopy() because it cuases problems.
	copy_thlen!(thlen1,thlen0)
	copy_thlen!(thlen2,thlen1)
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
	lenfun(len::Float64) = abs(len)^(-2)*cdfscale(abs(len))
	# The first value used in the Gaussian filter.
	sig1 = sqrt( lenfun(len1) + lenfun(len2) )
	sig1 *= 2*pi*sqrt(epsilon*dt)
	# The second value used in the Gaussian filter.
	sig2 = sqrt( lenfun(len0) + 2*lenfun(len1) + lenfun(len2) )
	sig2 *= 2*pi*sqrt(epsilon*dt)
	# Apply the appropriate Gaussian filters to advance theta in time.
	thlen2.theta = gaussfilter(theta1 - 2*pi*alpha, sig1) + 2*pi*alpha
	thlen2.theta += 0.5*dt*( 3*gaussfilter(n1,sig1) - gaussfilter(n0,sig2) )
	# Apply a Krasny-filter to theta (voided)
	### thlen2.theta = krasnyfilter(thlen2.theta - 2*pi*alpha) + 2*pi*alpha	
	return
end

#--------------- DETAILS OF THE MULTISTEP METHOD ---------------#
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
	

	# HERE
	nterm = (datau + dtheta.*vtang)/len

	# HERE?
	xsmdot = mean(-vnorm.*sin(theta) + vtang.*cos(theta))
	ysmdot = mean( vnorm.*cos(theta) + vtang.*sin(theta))
	return mterm, nterm, xsmdot, ysmdot
end
# tangvel: Compute the tangential velocity and mterm = dL/dt along the way.
function tangvel(dtheta::Vector{Float64}, vnorm::Vector{Float64})

	# HERE
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

#--------------- OTHER ---------------#
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
# getalpha: Calculate the parameterization variable, alpha = s/L, using an offset grid.
function getalpha(npts::Integer)
	dalpha = 1.0/npts
	return alpha = collect(range(0.5*dalpha, dalpha, npts))
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