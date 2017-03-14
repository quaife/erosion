# Test the curvature driven flow

function test()
	# Test convergence with respect to dt.
	npts = 256
	dt = 1e-2
	ndts = 3
	# Usually fixed parameters.
	epsilon = 0.5
	tfin = 0.1
	lenevo = 0
	len0 = 2*pi*0.2
	a1 = 0.1
	# Enter the time loop.
	L2ev = zeros(Float64,ndts)
	Linfev = zeros(Float64,ndts)
	thnew = cdf(npts,dt,epsilon,tfin,lenevo,len0,a1)
	for nn=1:ndts
		dt = 0.5*dt
		println("dt = ", dt)
		thold = thnew
		thnew = cdf(npts,dt,epsilon,tfin,lenevo,len0,a1)
		#L2ev[nn] = L2err(thold,thnew)
		Linfev[nn] = Linferr(thold,thnew)
	end
	#println(L2ev)
	println(Linfev)
end

# Curvature-driven flow.
function cdf(npts::Int, dt::Float64, epsilon::Float64,
		tfin::Float64, lenevo::Int, len0::Float64, a1::Float64)
	# Calculate parameters and initial shape.
	thlenvec0 = makeshape(npts,len0,a1)
	params = ParamType(dt,epsilon,0.,0.,lenevo,0)
	nsteps = round(Int,tfin/dt)
	# Start the simulation with RKstarter.
	thlenvec1 = RKstarter!(thlenvec0,params,noatau!)
	# Continue the simulation with the multi-step method.
	for nn = 2:nsteps
		noatau!(thlenvec1,params)
		advance_thetalen!(thlenvec1,thlenvec0,params)
	end
	thfin = thlenvec1[1].theta
	return thfin
end

# The Linf-difference between two functions
function Linferr(f1::Vector{Float64},f2::Vector{Float64})
	N1 = length(f1)
	N2 = length(f2)
	if N1 != N2
		throw("The vectors have to be equal length",N1,N2)
	else
		err = maximum(abs(f1-f2))/maximum(abs(f2))
	end
	return err

end

#= THESE ARE WRONG right now.
# The L2-difference between two theta-type functions.
function L2err_th(th1::Vector{Float64},th2::Vector{Float64})
	N1 = length(th1)
	N2 = length(th2)
	if N1>N2 
		throw("The second vector should be longer ",N1,N2)
	else

		err = L2norm([f1; zeros(N2-N1)] - f2)
	end
	return err
end
# The L2-difference between two periodic functions.
function L2err_per(f1::Vector{Float64},f2::Vector{Float64})
	N1 = length(f1)
	N2 = length(f2)
	if N1>N2 
		throw("The second vector should be longer ",N1,N2)
	else
		err = L2norm([f1; zeros(N2-N1)] - f2)
	end
	return err
end
=#

# L2-norm of a periodic function using Parseval's identity.
function L2norm(ff::Vector{Float64})
	fh = fft(ff)/length(ff)
	norm2 = sum(abs(fh).^2)
	norm2 = sqrt(norm2)
	return norm2
end

# Make the initial shape.
function makeshape(npts::Int, len::Float64, a1::Float64, 
		xsm::Float64=0., ysm::Float64=0.)
	alpha = getalpha(npts)
	# theta for a circle.
	theta = 0.5*pi + 2*pi*alpha[:]
	# Add a sin component.
	theta += a1*sin(2*pi*alpha[:])
	# Put the shape in a thlenvec.
	thlen = new_thlen()
	thlen.theta = theta
	thlen.len = len
	thlen.xsm = xsm
	thlen.ysm = ysm
	thlenvec = [thlen]
	return thlenvec
end

# Function to pass into RKstarter without computing abs tau.
function noatau!(thlenv::Vector{ThetaLenType}, params::ParamType)
	nbods = length(thlenv)
	npts = length(thlenv[1].theta)
	for nn=1:nbods
		thlenv[nn].atau = zeros(Float64,npts)
	end
	return
end

