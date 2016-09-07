# Erode multiple bodies.
include("basic.jl")
function main()
	##### PARAMETERS #####
	# Geometry parameters.
	npts = 256
	nbods = 4
	xsm1, ysm1 = +0.0, +0.4
	xsm2, ysm2 = -0.0, -0.4
	xsm3, ysm3 = +0.4, +0.0
	xsm4, ysm4 = -0.4, -0.0	
	# For circle geometry
	rad1, rad2, rad3, rad4 = 0.2, 0.2, 0.2, 0.2
	# Evolution parameters.
	dt = 5e-4
	epsilon = 5e-2
	sigma = epsilon
	tfin = 10*dt
	# Misc parameters.
	beta = 0
	axlim = 1.0
	######################

	# Put the parameters in a single variable.
	params = ParamType(dt,epsilon,sigma,beta)
	# Create the initial geometries.
	thlen01 = circgeo(npts,rad1,xsm1,ysm1)
	thlen02 = circgeo(npts,rad2,xsm2,ysm2)
	thlen03 = circgeo(npts,rad3,xsm3,ysm3)
	thlen04 = circgeo(npts,rad4,xsm4,ysm4)
	# Create the vector of ThetaLenType values.
	thlenvec0 = [thlen01, thlen02, thlen03, thlen04]
	#thlenvec0 = [thlen01]
	
	# Plot the initial geometries.
	plotcurves!(thlenvec0, 0, axlim=axlim)	
	# Use RK2 as a starter.
	thlenvec1 = RKstarter!(thlenvec0, params)
	
	# Initialize values for the while loop (with slight adjustment to tfin).
	tfin -= 0.5*dt; tm = dt; cnt = 1; 
	# Plot the result for t=dt.
	plotcurves!(thlenvec1, cnt, axlim=axlim)
	# Enter while loop.
	while(tm < tfin)
		# Compute the new stress and save it.
		stokes!(thlenvec1, sigma)
		# Advance thlen forward in time using the multi-step method.
		advance_thetalen!(thlenvec1,thlenvec0,params)
		# Advance time & counter and plot the result.
		tm += dt; cnt += 1; plotcurves!(thlenvec1,cnt,axlim=axlim)
	end
end

main()
