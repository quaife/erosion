# Erode a body.
function main()
	include("basic.jl")
	##### PARAMETERS #####
	# Geometry parameters.
	npts = 256
	xsm = 0.5; ysm = 0.4;
	rad = 0.3								# For circle geometry.
	nsides = 4; sigma = 0.1; sdlen = 0.2	# For polygon geometry.
	# Evolution parameters.
	tfin = 0.1
	dt = 2e-4
	epsilon = 0.2
	beta = 0
	# Misc parameters.
	axlim = 1.0
	######################

	# Make slight adjustment to ensure that tfin is obtained.
	tfin += 0.5*dt
	# Put the parameters in a single variable.
	params = ParamType(dt,epsilon,beta)
	# Create the initial geometry.
	thlen00 = circgeo(npts,rad,xsm,ysm)
	#thlen00 = polygongeo(npts, nsides,sigma,sdlen,xsm,ysm)
	# Copy to thlen0, which will be modified in the multi-step method.
	thlen0 = new_thlen()
	copy_thlen!(thlen00,thlen0)
	# Use RK2 as a starter.
	thlen1 = RKstarter!(thlen0,params)
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
