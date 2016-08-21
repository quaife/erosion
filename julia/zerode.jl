# Erode a body.
function main()
	include("basic.jl")
	##### PARAMETERS #####
	# Geometry parameters.
	npts = 256
	rad = 0.2
	# Evolution parameters.
	tfin = 0.1
	dt = 1e-3
	epsilon = 0.02
	beta = 0
	######################

	# Make slight adjustment to ensure that tfin is obtained.
	tfin += 0.5*dt
	# Put the parameters in a single variable.
	params = ParamType(dt,epsilon,beta)
	# Create the initial circular geometry.
	thlen00 = circgeo(npts,rad)
	# Copy to thlen0, which will be modified in the multi-step method.
	thlen0 = new_thlen()
	copy_thlen!(thlen00,thlen0)
	# Use RK2 as a starter.
	thlen1 = RKstarter!(thlen0,params)
	# Plot the result.
	tm = dt; cnt = 1; plotcurve!(thlen1,thlen00,cnt)
	# Enter while loop.
	while(tm < tfin)
		# Compute the new stress and save it.
		stokes!([thlen1])
		# Advance thlen forward in time using the multi-step method.
		advance_thetalen!(thlen1,thlen0,params)
		# Advance time & counter and plot the result.
		tm += dt; cnt += 1; plotcurve!(thlen1,thlen00,cnt)
	end
end

main()
