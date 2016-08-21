#= misc.jl: Miscellaneous functions =#

# ParamType includes the parameters dt, epsilon, and beta.
type ParamType
	dt::Float64; epsilon::Float64; beta::Real;
end
# ThetaLenType includes all of the data that for a curve.
type ThetaLenType
	theta::Vector{Float64}; len::Float64; xc::Float64; yc::Float64; 
	atau::Vector{Float64}; mterm::Float64; nterm::Vector{Float64};
	xx::Vector{Float64}; yy::Vector{Float64};
end
# Create a new ThetaLenType that has all zeros.
function new_thlen()
	return ThetaLenType([0.0], 0.0, 0.0, 0.0, [0.0], 0.0, [0.0], [0.0], [0.0])
end
# Copy the relevant contents from thlen1 to thlen2.
function copy_thlen!(thlen1::ThetaLenType, thlen2::ThetaLenType)
	thlen2.theta = thlen1.theta
	thlen2.len = thlen1.len
	thlen2.xc = thlen1.xc
	thlen2.yc = thlen1.yc
	thlen2.atau = thlen1.atau
	thlen2.mterm = thlen1.mterm
	thlen2.nterm = thlen1.nterm
	thlen2.xx = thlen1.xx
	thlen2.yy = thlen1.yy
end


########## Starter routines ##########
# RKstarter: Explicit second-order Runge-Kutta to start the time stepping.
function RKstarter!(thlen0::ThetaLenType, params::ParamType)
	# Extract the needed variables.
	dt, epsilon, beta = params.dt, params.epsilon, params.beta
	theta0, len0 = thlen0.theta, thlen0.len
	# Get the time derivatives at t=0.
	thlen0.atau = stokes_thl_sing(theta0,len0)
	th0dot, m0 = thetadot!(thlen0,params)
	# Take the first half-step of RK2.
	len05 = len0 + 0.5*dt*m0
	theta05 = theta0 + 0.5*dt*th0dot
	# Get the time derivatives at t=0.5*dt (do not need to use ThetaLenType).
	atau05 = stokes_thl_sing(theta05,len05)
	th05dot,m05,n05 = thetadot(theta05,len05,atau05,params)
	# Create a new ThetaLenType variables and take the second step of RK2.
	thlen1 = new_thlen()
	thlen1.len = len0 + dt*m05
	thlen1.theta = theta0 + dt*th05dot
	return thlen1
end

# thetadot: Calculate the time derivative of theta.
function thetadot(theta::Vector{Float64}, len::Float64, atau::Vectlor{Float64}, params::ParamType)
	# Extract the needed variables.
	dt, epsilon, beta = params.dt, params.epsilon, params.beta
	alpha = getalpha(endof(theta))
	# Calculate mterm and nterm at time 0.
	mterm, nterm = getmn(theta,len,atau,params)
	# Calculate the time derivative of theta.
	dth = specdiff(theta - 2*pi*alpha) + 2*pi
	d2th = specdiff(dth)
	thdot = epsilon*len^(beta-2)*d2th + nterm
	# Return thdot, mterm, nterm.
	return thdot, mterm, nterm
end
#= thetadot: Dispatch for ThetaLenType.
It also calculates mterm and nterm and saves them in thlen. =#
function thetadot!(thlen::ThetaLenType, params::ParamType)
	# Call thetadot and save mterm and nterm in thlen.
	thdot, thlen.mterm, thlen.nterm = thetadot(thlen.theta, thlen.len, thlen.atau, params)
	return thdot, mterm
end
########################################






#= getxy: Given theta and len, reconstruct the x and y coordinates of a body.
xc and yc are the coordinates of the center of mass (???).
While we're at it, also calculate the normal direcations. =#
function getxy(theta::Vector{Float64}, len::Float64, xc::Float64, yc::Float64)
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
# Create another dispatch of getxy for input of type ThetaLenType.
function getxy!(thlen::ThetaLenType)
	thlen.xx, thlen.yy = getxy(thlen.theta, thlen.len, thlen.xc, thlen.yc)
end

# plotcurve: Plot a curve from the theta-len values.
function plotcurve(thlen::ThetaLenType, x0::Vector{Float64}, y0::Vector{Float64}, cnt::Integer)
	# Set the axis limits.
	axlim = 0.5
	# Reconstruct the x,y coordinates of the curve.
	xx,yy = getxy(thlen)
	# Make the plot.
	p1 = plot(x0,y0,"-", xx,yy,"-")
	xlim(-axlim,axlim); ylim(-axlim,axlim)
	# Save the figures in a folder.
	figname = string("../figs/fig",string(cnt),".pdf")
	savefig(p1, figname, width=500, height=500)
end
