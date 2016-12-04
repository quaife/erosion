# main.jl

# erosion: The main routine to erode a group of bodies for input of Vector{ThetaLenType}.
function erosion(tfin::Float64, dt::Float64, thlenvec0::Vector{ThetaLenType}; 
	axlims::Vector{Float64} = [1.,1.])
	# Extract the basic parameters
	npts = endof(thlenvec0[1].theta)
	nbods = endof(thlenvec0)
	nsteps = round(Int,tfin/dt)
	# Calculate the smoothing parameters based on the spatial resolution.
	epsilon = 10./npts
	sigma = 15./npts
	params = ParamType(dt,epsilon,sigma,0)
	# Set up the target points to measure u,v,p.
	ntar0 = 10; xmax = 2.8; ymax = 0.8
	ntar,xtar,ytar,utar,vtar,ptar = targets(ntar0,xmax,ymax)

	# Use the Runge-Kutta starter, while plotting and saving before and after.
	plotnsave(thlenvec0,0,axlims=axlims)
	thlenvec1 = RKstarter!(thlenvec0, params)
	plotnsave(thlenvec1,1,axlims=axlims)
	# Enter the time loop.
	for cnt = 2:nsteps
		utar,vtar,ptar = stokes!(thlenvec1,sigma,ntar,xtar,ytar)
		advance_thetalen!(thlenvec1,thlenvec0,params)
		plotnsave(thlenvec1,cnt,axlims=axlims)
	end
	return
end
# targets: Set up the target points to measure velocity and pressure: u,v,p.
function targets(nn::Integer, xmax::Float64, ymax::Float64)
	# Make the grid.
	ytar = collect(linspace(-ymax,ymax,nn))
	ytar = [ytar; ytar]
	xtar = ones(Float64,nn)
	xtar = xmax*[-xtar; xtar]
	# Initialize u,v,p at target points.
	utar,vtar,ptar = [zeros(Float64,2*nn) for ii=1:3]
	return 2*nn,xtar,ytar,utar,vtar,ptar
end

#################### Starter routines ####################
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
# festep: Dispatch for ThetaLenType.
# Allow the starting point and the point where the derivatives are taken to be different.
function festep(dt::Float64, thdot::Vector{Float64}, 
		thlen0::ThetaLenType, thlendots::ThetaLenType)
	thlen1 = new_thlen()
	thlen1.theta, thlen1.len, thlen1.xsm, thlen1.ysm = festep(dt, thdot, 
		thlen0.theta, thlen0.len, thlen0.xsm, thlen0.ysm, 
		thlendots.mterm, thlendots.xsmdot, thlendots.ysmdot)
	return thlen1
end
#= RKstarter!: Explicit second-order Runge-Kutta to start the time stepping.
Works for vectors of ThetaLenType.
It also calculates mterm, nterm, xsmdot, ysmdot and saves them in thlenvec0. =#
function RKstarter!(thlenvec0::Vector{ThetaLenType}, params::ParamType)
	dt = params.dt
	sigma = params.sigma
	nbods = endof(thlenvec0) 
	thlenvec05 = [new_thlen() for ii=1:nbods]
	thlenvec1 = [new_thlen() for ii=1:nbods]
	# Compute the stress at t=0 and take the first step of RK2
	stokes!(thlenvec0, sigma)
	for ii = 1:nbods
		thlen0 = thlenvec0[ii]
		thdot = thetadot!(thlen0,params)
		thlenvec05[ii] = festep(0.5*dt, thdot, thlen0, thlen0)	
	end
	# Compute the stress at t=0.5*dt and take the second step of RK2
	stokes!(thlenvec05, sigma)	
	for ii = 1:nbods
		thlen0 = thlenvec0[ii]
		thlen05 = thlenvec05[ii]
		thdot = thetadot!(thlen05,params)
		thlenvec1[ii] = festep(dt, thdot, thlen0, thlen05)		
	end
	# Remove curves with non-positive length and return.
	trimthlenvec!(thlenvec1, thlenvec0)
	return thlenvec1
end
