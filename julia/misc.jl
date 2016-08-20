#= misc.jl: Miscellaneous functions =#

########## Starter routines ##########
# RKstarter: Explicit second-order Runge-Kutta to start the time stepping.
function RKstarter(theta0::Vector{Float64}, len0::Float64,
		dt::Float64, epsilon::Float64, beta::Real)
	# Get the time derivatives at t=0.
	atau = stokes_thl_sing(theta0,len0)
	th0dot,M0,N0 = thetadot(atau,theta0,len0,epsilon,beta)
	# Take the first half-step of RK2.
	len05 = len0 + 0.5*dt*M0
	theta05 = theta0 + 0.5*dt*th0dot
	# Get the time derivatives at t=0.5*dt.
	atau = stokes_thl_sing(theta05,len05)
	th05dot,M05,N05 = thetadot(atau,theta05,len05,epsilon,beta)
	# Take the second step of RK2.
	len1 = len0 + dt*M05
	theta1 = theta0 + dt*th05dot
	return theta1,len1,M0,N0
end
#= thetadot: Calculate the time derivative of theta;
Also return MM and NN while we're at it. Only used in the RKstarter. =#
function thetadot(atau::Vector{Float64}, theta::Vector{Float64}, len::Float64, 
		epsilon::Float64, beta::Real)
	# Get the M and N terms
	MM,NN = getMN(atau,theta,len,epsilon,beta)
	# Calculate the time derivative of theta.
	alpha = getalpha(endof(theta))
	dth = specdiff(theta - 2*pi*alpha) + 2*pi
	d2th = specdiff(dth)
	thdot = epsilon*len^(beta-2)*d2th + NN
	return thdot, MM, NN
end
########################################

#= getxy: Given theta and len, reconstruct the x and y coordinates of a body.
xc and yc are the coordinates of the center of mass (???).
While we're at it, also calculate the normal direcations. =#
function getxy(theta::Vector{Float64}, len::Float64, xc::Float64=0.0, yc::Float64=0.0)
	# The increments of dx and dy
	dx = len * (cos(theta) - mean(cos(theta)))
	dy = len * (sin(theta) - mean(sin(theta)))
	# Integrate to get the x,y coordinates.
	xx = specint(dx)
	yy = specint(dy)

	## TO DO: Move the center of mass or average values.

	# The normal vector, direction???
	nx = -sin(theta)
	ny = cos(theta)
	return xx,yy
end

# plotcurve: Plot a curve from the theta-len values.
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
