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

#=

# Enter while loop.
while(tm < tfin)
	# Compute the stress
	atau = stokes_thl_sing(theta1,len1)
	# Advance theta and len in time
	theta1,len0,len1,M0,N0 = advance_thetalen(atau,theta1,
								len0,len1,M0,N0,params)
	# Advance time and counter, and plot the result.
	tm += dt; cnt += 1; plotcurve(theta1,len1,x0,y0,cnt)
end
=#