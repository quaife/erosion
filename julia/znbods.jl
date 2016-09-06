# Erode a body.
include("basic.jl")

#= RKstarter: Explicit second-order Runge-Kutta to start the time stepping.
It also calculates mterm and nterm and saves them in thlen0. =#
function RKstarter!(thlenvec0::Vector{ThetaLenType}, params::ParamType)
	# Extract the needed variables.
	dt, epsilon, beta = params.dt, params.epsilon, params.beta
	# Get the time derivatives at t=0.
	stokes!(thlenvec0)

	th0dot, m0 = thetadot!(thlenvec0,params)



	# Create a new ThetaLenType variables for t=0.5*dt.
	thlen05 = new_thlen()
	# Take the first half-step of RK2.
	thlen05.len = len0 + 0.5*dt*m0
	thlen05.theta = theta0 + 0.5*dt*th0dot
	# Update the surface-mean coordinates.
	thlen05.xsm = thlen0.xsm + 0.5*dt*thlen0.xsmdot
	thlen05.ysm = thlen0.ysm + 0.5*dt*thlen0.ysmdot
	# Get the time derivatives at t=0.5*dt.
	stokes!([thlen05])
	th05dot, m05 = thetadot!(thlen0,params)

	# Create a new ThetaLenType variables for time t=dt.
	thlen1 = new_thlen()
	# Take the second step of RK2.
	thlen1.len = len0 + dt*m05
	thlen1.theta = theta0 + dt*th05dot
	# Update the surface-mean coordinates.
	thlen1.xsm = thlen0.xsm + dt*thlen05.xsmdot
	thlen1.ysm = thlen0.ysm + dt*thlen05.ysmdot
	# Return
	return thlen1
end

# Main program.
function main()
	##### PARAMETERS #####
	# Geometry parameters.
	npts = 256
	nbods = 2
	xsm1, ysm1 =  0.4, 0.0
	xsm2, ysm2 = -0.4, 0.0
	# For circle geometry
	rad1, rad2 = 0.2, 0.2
	# Evolution parameters.
	tfin = 1e-2
	dt = 5e-4
	epsilon = 0.2
	beta = 0
	# Misc parameters.
	axlim = 1.0
	######################

	# Create the initial geometries.
	thlen001 = circgeo(npts,rad1,xsm1,ysm1)
	thlen002 = circgeo(npts,rad2,xsm2,ysm2)
	# Create the vector of ThetaLenType values.
	thlenvec0 = [thlen001, thlen002]

	# Make slight adjustment to ensure that tfin is obtained.
	tfin += 0.5*dt
	# Put the parameters in a single variable.
	params = ParamType(npts,nbods,dt,epsilon,beta)
	# Use RK2 as a starter.
	thlenvec1 = RKstarter!(thlenvec0,params)




	# Plot the result.
	tm = dt; cnt = 1; plotcurve!(thlen1,thlen00,cnt,axlim=axlim)
	# Enter while loop.
	while(tm < tfin)
		# Compute the new stress and save it.
		stokes!([thlen1])
		# Advance thlen forward in time using the multi-step method.
		advance_thetalen!(thlen1,thlen0,params)
		# Advance time & counter and plot the result.
		tm += dt; cnt += 1; plotcurve!(thlen1,thlen00,cnt,axlim=axlim)
	end
end

main()
