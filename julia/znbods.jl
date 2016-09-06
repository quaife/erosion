# Erode multiple bodies.
include("basic.jl")
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

	# Put the parameters in a single variable.
	params = ParamType(dt,epsilon,beta)
	# Create the initial geometries.
	thlen001 = circgeo(npts,rad1,xsm1,ysm1)
	thlen002 = circgeo(npts,rad2,xsm2,ysm2)
	# Create the vector of ThetaLenType values.
	thlenvec0 = [thlen001, thlen002]
	# Plot the initial geometries.
	plotcurves!(thlenvec0,0,axlim=axlim)	
	# Use RK2 as a starter.
	thlenvec1 = RKstarter!(thlenvec0,params)
	
	# Initialize values for the while loop (with slight adjustment to tfin).
	tfin += 0.5*dt; tm = dt; cnt = 1; 
	# Plot the result for t=dt.
	plotcurves!(thlenvec1,cnt,axlim=axlim)
	# Enter while loop.
	while(tm < tfin)
		# Compute the new stress and save it.
		stokes!(thlenvec1)
		# Advance thlen forward in time using the multi-step method.
		advance_thetalen!(thlenvec1,thlenvec0,params)
		# Advance time & counter and plot the result.
		tm += dt; cnt += 1; plotcurves!(thlenvec1,cnt,axlim=axlim)
	end
end

main()