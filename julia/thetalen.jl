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



function RK2(thld0::ThLenDenType)
	getstress!(thld0, params)
	thld05 = feuler(thld0, 0.5*dt, thld0)
	getstress!(thld05, params)
	thld1 = feuler(thld0, dt, thld05)
end

function feuler(thld0::ThLenDenType, dt, derivs::ThLenDenType)
	npts, nbods = getnvals(thld0.thlenvec)
	thlv1 = new_thlenvec(npts)
	for nn = 1:nbods
		thlv1[nn] = feuler(thld0.thlenvec[nn], dt, derivs.thlenvec[nn].atau)
	end
	thld1 = new_thlenden(thlv1)
	return thld1
end

function feuler(thl0::ThetaLenType, dt::Float64, atau::Vector{Float64})
	mterm, nterm, xsmdot, ysmdot = getmn( ??? )
	len05 = thl0.len + dt*mterm

	theta05 = gaussfilter(thl0.theta + dt*nterm)













#--------------- SMALL ROUTINES ---------------#
# getmn: Calculates mterm and nterm: mterm=dL/dt and nterm is the nonlinear term.
function getmn(theta::Vector{Float64}, len::Float64, atau::Vector{Float64}, epsilon::Float64)
	if atau==[]; throw("atau has not been computed"); return; end
	alpha = getalpha(endof(theta))
	dtheta = specdiff(theta - 2*pi*alpha) + 2*pi
	vnorm = atau + epsilon*cdfscale(len)*len^(-1) * (dtheta - 2*pi)
	vtang, mterm = tangvel(dtheta, vnorm)
	datau = specdiff(atau)	# Derivative of absolute-value of shear stress.
	nterm = (datau + vecmult(dtheta,vtang))/len
	# Get the derivatives of xsm and ysm: should I dealias products here?
	xsmdot = mean(-vnorm.*sin(theta) + vtang.*cos(theta))
	ysmdot = mean( vnorm.*cos(theta) + vtang.*sin(theta))
	return mterm, nterm, xsmdot, ysmdot
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
