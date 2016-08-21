# Erode a body.
include("includes.jl")

##### PARAMETERS #####
# Geometry parameters.
npts = 128
rad = 0.2
# Evolution parameters.
tfin = 0.5
dt = 0.001
epsilon = 0.02
beta = 0
######################

# Make slight adjustment to ensure that tfin is obtained.
tfin += 0.5*dt
# Put the parameters in a single variable.
params = ParamType(dt,epsilon,beta)
# Create the initial circular geometry.
thlen0 = circgeo(npts,rad)
# Get the initial x and y coordinates
x0,y0 = getxy(thlen0)
# Use RK2 as a starter.
thlen1 = RKstarter(thlen0,params)

# Plot the result.
tm = dt; cnt = 1; plotcurve(theta1,len1,x0,y0,cnt)


# Enter while loop.
while(tm < tfin)
	# Compute the new stress and save it.
	atau = stokes_thl_sing(theta1,len1)
	thlen1.atau = atau
		# Advance thlen forward in time using the multi-step method.
	advance_thetalen!(thlen1,thlen0,params)
	# Advance time & counter and plot the result.
	tm += dt; cnt += 1; plotcurve(theta1,len1,x0,y0,cnt)
end
