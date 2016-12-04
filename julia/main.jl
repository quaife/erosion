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
	# Save the parameters in one variable.
	params = ParamType(dt,epsilon,sigma,0)
	# The output file to save the geometry.
	outfile = "geoout.dat"
	iostream = open(outfile,"w")
	writedlm(iostream,[npts,nbods,nsteps+1])
	close(iostream)


	# TEMPORARY: Set up the target points to measure u, v, and pressure.
	ntargs = 11
	yend = 0.8
	ytar = collect(linspace(-yend, yend, ntargs))
	xtar = -2.8 * ones(Float64, ntargs)
	# Initialize variables for u, v, and pressure at target points.
	utar,vtar,ptar = [zeros(Float64,ntargs) for ii=1:3]
	pavg = zeros(Float64,nsteps)


	# Plot the initial geometries, t=0, and save in a data file.
	plotcurves!(thlenvec0,0; axlims=axlims)
	savexydata(outfile, thlenvec0)
	# Use RK2 as a starter.
	thlenvec1 = RKstarter!(thlenvec0, params)
	# Plot the result for t=dt.
	plotcurves!(thlenvec1,1; axlims=axlims)
	savexydata(outfile, thlenvec1)
	# Enter the time loop.
	for cnt = 2:nsteps
		# Compute the new stress and save it.
		utar,vtar,ptar = stokes!(thlenvec1, params.sigma, ntargs, xtar, ytar)
		# Advance thlen forward in time using the multi-step method.
		advance_thetalen!(thlenvec1,thlenvec0,params)
		# Calculate the average pressure.
		pavg[cnt] = mean(ptar)
		# Plot the results.
		plotcurves!(thlenvec1,cnt; axlims=axlims)
		#plotpress(ytar,ptar,cnt)

		# Save the current geometry in a data file.
		savexydata(outfile, thlenvec1)
	end
	return
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
	# Initialize vectors of ThetaLenType.
	thlenvec05 = [new_thlen() for ii=1:nbods]
	thlenvec1 = [new_thlen() for ii=1:nbods]
	# Compute the stress at t=0.
	stokes!(thlenvec0, sigma)
	# For each body, take the first step of RK2.
	for ii = 1:nbods
		# Need thlen0 for each body.
		thlen0 = thlenvec0[ii]
		# Calculate the time derivatives: thdot, mterm, xsmdot, ysmdot.
		thdot = thetadot!(thlen0,params)
		# Take the first step of RK2.
		thlenvec05[ii] = festep(0.5*dt, thdot, thlen0, thlen0)	
	end
	# Compute the stress at t=0.5*dt.
	stokes!(thlenvec05, sigma)	
	# For each body, take the second step of RK2.
	for ii = 1:nbods
		# Need both thlen0 and thlen05 for each body.
		thlen0 = thlenvec0[ii]
		thlen05 = thlenvec05[ii]
		# Calculate the time derivatives: thdot, mterm, xsmdot, ysmdot.
		thdot = thetadot!(thlen05,params)
		# Take the second step of RK2.
		thlenvec1[ii] = festep(dt, thdot, thlen0, thlen05)		
	end
	return thlenvec1
end
