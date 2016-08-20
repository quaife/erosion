# Test curvature-driven flow for an initial triangle.
include("includes.jl")

# plotcurve
function plotcurve(theta::Vector{Float64}, len::Float64, 
		x0::Vector{Float64}, y0::Vector{Float64}, cnt::Integer)
	# Reconstruct the x,y coordinates of the curve.
	xx,yy = getxy(theta,len)
	p1 = plot(x0,y0,".-", xx,yy,".-")
	xlim(-1.0,1.0); ylim(-1.0,1.0)
	figname = string("../figs/fig",string(cnt),".pdf")
	savefig(p1, figname, width=500, height=500)
end

##### PARAMETERS #####
# Triangle geometry parameters.
nn = 128
angle = 90
sigma = 0.05
# Evolution parameters.
tfin = 2.0
dt = 0.01
epsilon = 0.1
beta = 0
######################

# Get the initial triangular geometry.
theta0,len0,xback = trigeo(nn,angle,sigma)
x0,y0 = getxy(theta0,len0)
# For pure curvature driven flow, set the stress term to zero.
atau = zeros(Float64,nn)

# Initialize with RK2.
# Get the time derivatives at t=0.
th0dot,M0,N0 = thetadot(atau,theta0,len0,epsilon,beta)
# Take the first half-step of RK2.
len05 = len0 + 0.5*dt*M0
theta05 = theta0 + 0.5*dt*th0dot
# Get the time derivatives at t=0.5*dt.
th05dot,M05,N05 = thetadot(atau,theta05,len05,epsilon,beta)
# Take the second step of RK2.
len1 = len0 + dt*M05
theta1 = theta0 + dt*th05dot
# Plot the result.
cnt = 1; plotcurve(theta1,len1,x0,y0,cnt)

# Enter while loop.
tfin += 0.1*dt
tm = dt
while(tm < tfin)
	# Advance theta and len in time
	theta1,len0,len1,M0,N0 = advance_thetalen(atau,theta1,
								len0,len1,M0,N0,dt,epsilon,beta)
	# Advance time and counter, and plot the result.
	tm += dt; cnt += 1; plotcurve(theta1,len1,x0,y0,cnt)
end

# Plot area vs. time