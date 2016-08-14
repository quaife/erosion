# Test curvature-driven flow for an initial triangle.
# workspace()
include("spectral.jl")
include("thetalen.jl")
include("geometries.jl")
using Winston

# plotcurve
function plotcurve(theta::Vector{Float64}, len::Float64, 
		xback::Float64, x0::Vector{Float64}, y0::Vector{Float64}, cnt::Integer)
	# Reconstruct the x,y coordinates of the curve.
	xx,yy = getcurve(theta,len,xback)
	p1 = plot(x0,y0,".-", xx,yy,".-")
	xlim(-1.0,1.0); ylim(-1.0,1.0)
	figname = string("../Figs/fig",string(cnt),".pdf")
	savefig(p1, figname, width=500, height=500)
end

##### PARAMETERS #####
# Triangle geometry parameters.
nn = 128
angle = 90
sigma = 0.05
# Evolution parameters.
tfin = 2.0
dt = 0.02
epsilon = 0.1
beta = 0
######################

# Get the initial triangular geometry.
theta,len,xback,alpha = trigeo(nn, angle, sigma)
x0,y0 = getcurve(theta,len,xback)
# For pure curvature driven flow, make the stress term vanish.
ast = 0.0*theta[:]

# Initialize with RK2
len0 = len
# Get the time derivatives at t=0
M0,N0 = getMN(ast,theta,len,epsilon,beta)
th0dot = thetadot(theta,len,N0,epsilon,beta)
# Take the first half step of RK2
len05 = len0 + 0.5*dt*M0
th05 = theta + 0.5*dt*th0dot
# Get the time derivatives at t=0.5*dt
M05,N05 = getMN(ast,th05,len05,epsilon,beta)
th05dot = thetadot(th05,len05,N05,epsilon,beta)
# Take the second half step of RK2
len1 = len05 + 0.5*dt*M05
theta1 = th05 + 0.5*dt*th05dot
# Plot the result.
cnt = 1
plotcurve(theta1,len1,xback,x0,y0,cnt)

# Enter while loop.
tfin += 0.1*dt
tm = dt
while(tm < tfin)
	# Advance theta and len in time
	theta1,len0,len1,M0,N0 = advance_thetalen(ast,theta1,
								len0,len1,M0,N0,dt,epsilon,beta)
	# Advance time and counter.
	tm += dt
	cnt += 1
	# Plot the result.
	plotcurve(theta1,len1,xback,x0,y0,cnt)
end
