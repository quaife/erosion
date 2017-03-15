# Test the curvature driven flow

function test()
	# Test convergence with respect to dt.
	npts = 256
	dt = 5e-3
	ndts = 3
	# Usually fixed parameters.
	epsilon = 0.1
	tfin = 0.2
	lenevo = 0
	len0 = 2*pi*0.2
	a1 = 0.5
	
	# Initialization.
	thlenvec = makeshape(npts,len0,a1)
	# Make plot of initial shape.
	plotfolder = "../figs/"
	newfolder(plotfolder)
	plotshapetheta(thlenvec,plotfolder,0)
	# Run the first simulation.
	thlenvec = cdf(npts,dt,epsilon,tfin,lenevo,len0,a1)
	thnew = thlenvec[1].theta
	# Enter the time loop to run many simulations and calculate errors.
	L2v = zeros(Float64,ndts)
	Linfv = zeros(Float64,ndts)
	for nn=1:ndts
		dt = 0.5*dt
		thold = thnew
		# Run the simulation.
		thlenvec = cdf(npts,dt,epsilon,tfin,lenevo,len0,a1)
		thnew = thlenvec[1].theta
		# Calculate errors.
		Linfv[nn] = Linferr(thold,thnew)
		L2v[nn] = L2err_th(thold,thnew)
		plotshapetheta(thlenvec,plotfolder,nn)
	end
	orderinf = order(2,Linfv)
	order2 = order(2,L2v)
	println(Linfv)
	println(L2v)
	println(orderinf)
	println(order2)
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
	return thlenvec1
end

function plotshapetheta(thlenvec::Vector{ThetaLenType}, plotfolder::AbstractString, nn::Int)
	figshape = string(plotfolder,"shape",nn,".pdf")
	figtheta = string(plotfolder,"theta",nn,".pdf")
	plotcurves(thlenvec, figshape)
	# Plot theta.
	height = 400
	width = 600
	pp = plot()
	xlim(0,1.); ylim(0,3*pi)
	theta = thlenvec[1].theta
	npts = length(theta)
	alpha = getalpha(npts)
	pp = oplot(alpha,theta,"-")
	savefig(pp, figtheta, width=width, height=height)
end

#= Calculate the order of convergence given a vector of errors. =#
function order(Nfac::Integer, errv::Vector{Float64})
	return log(errv[1:end-1]./errv[2:end])./log(Nfac)
end
# The Linf-difference between two functions
function Linferr(f1::Vector{Float64},f2::Vector{Float64})
	if length(f1) != length(f2)
		throw("The vectors should be equal length")
	else
		err = maximum(abs(f1-f2))/maximum(abs(f2))
	end
	return err
end
# The L2-difference between two theta-type functions.
function L2err_th(th1::Vector{Float64},th2::Vector{Float64})
	# Make th1 and th2 periodic by removing the linear part.
	thp1 = th1 - 2*pi*getalpha(length(th1))
	thp2 = th2 - 2*pi*getalpha(length(th2))
	# Calculate the error of the periodic components.
	err = L2err_per(thp1,thp2)
	return err
end
# The L2-difference between two periodic functions.
function L2err_per(f1::Vector{Float64},f2::Vector{Float64})
	N1 = length(f1)
	N2 = length(f2)
	if N1>N2 
		throw("The second vector should be longer ",N1,N2)
	else
		fh1 = fft(f1)/N1
		fh2 = fft(f2)/N2
		err = L2norm(fh1-fh2)
	end
	return err
end
# L2-norm of a periodic function using Parseval's identity.
function L2norm{T<:Number}(fh::Vector{T})
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
	# Add a perturbation.
	theta += a1*sin(4*pi*alpha[:])
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

test()

