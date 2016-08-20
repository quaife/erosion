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



# plotcurve
function plotcurve(theta::Vector{Float64}, len::Float64, 
		x0::Vector{Float64}, y0::Vector{Float64}, cnt::Integer)
	# Set the axis limits.
	axlim = 0.5
	# Reconstruct the x,y coordinates of the curve.
	xx,yy = getxy(theta,len)
	# Make the plot.
	p1 = plot(x0,y0,"-", xx,yy,"-")
	xlim(-axlim,axlim); ylim(-axlim,axlim)
	# Save the figures in a folder.
	figname = string("../figs/fig",string(cnt),".pdf")
	savefig(p1, figname, width=500, height=500)
end

########## MAIN ##########
# Make slight adjustment to ensure that tfin is obtained.
tfin += 0.1*dt
# Create the initial circular geometry.
theta0,len0 = circgeo(npts,rad)
x0,y0 = getxy(theta0,len0)
# Use RK2 as a starter.
theta1,len1,M0,N0 = RKstarter(theta0,len0,dt,epsilon,beta)
# Plot the result.
tm = dt; cnt = 1; plotcurve(theta1,len1,x0,y0,cnt)

# Enter while loop.
while(tm < tfin)
	# Compute the stress
	atau = stokes_thl_sing(theta1,len1)
	# Advance theta and len in time
	theta1,len0,len1,M0,N0 = advance_thetalen(atau,theta1,
								len0,len1,M0,N0,dt,epsilon,beta)
	# Advance time and counter, and plot the result.
	tm += dt; cnt += 1; plotcurve(theta1,len1,x0,y0,cnt)
end
