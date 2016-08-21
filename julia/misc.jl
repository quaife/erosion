#= misc.jl: Miscellaneous functions =#

# ParamType includes the parameters dt, epsilon, and beta.
type ParamType
	dt::Float64; epsilon::Float64; beta::Real
end
# ThetaLenType includes all of the data that for a curve.
type ThetaLenType
	npts::Integer; alpha::Vector{Float64}; theta::Vector{Float64}; len::Float64
	xc::Float64; yc::Float64; atau::Vector{Float64}; mterm::Float64; nterm::Vector{Float64}
end
# Create a new ThetaLenType variables that only inherits npts and alpha.
function new_thlen(thlen::ThetaLenType)
	npts = thlen.npts
	zvec = zeros(Float64,npts)
	return ThetaLenType(npts, thlen.alpha, zvec, 0.0, 0.0, 0.0, zvec, 0.0, zvec)
end
# Copy the relevant contents from thlen1 to thlen2.
function copy_thlen!(thlen1::ThetaLenType, thlen2::ThetaLenType)
	thlen2.npts = thlen1.npts
	thlen2.alpha = thlen1.alpha
	thlen2.theta = thlen1.theta
	thlen2.len = thlen1.len
	thlen2.xc = thlen1.xc
	thlen2.yc = thlen1.yc
	thlen2.atau = thlen1.atau
	thlen2.mterm = thlen1.mterm
	thlen2.nterm = thlen1.nterm
end


function foo!(thlen::ThetaLenType)
	thlen.len = 55.0
	return
end


########## Starter routines ##########
# RKstarter: Explicit second-order Runge-Kutta to start the time stepping.
function RKstarter!(thlen0::ThetaLenType, params::ParamType)
	# Extract the needed variables.
	dt, epsilon, beta = params.dt, params.epsilon, params.beta
	theta0, len0 = thlen0.theta, thlen0.len
	# Get the time derivatives at t=0.
	thlen0.atau = stokes_thl_sing(theta0,len0)
	th0dot = thetadot!(thlen0,params)
	m0 = thlen0.mterm
	# Take the first half-step of RK2.
	len05 = len0 + 0.5*dt*m0
	theta05 = theta0 + 0.5*dt*th0dot
	# Save the results in a new ThetaLenType variable.
	thlen05 = new_thlen(thlen0)
	thlen05.theta = theta05
	thlen05.len = len05
	# Get the time derivatives at t=0.5*dt.
	thlen05.atau = stokes_thl_sing(theta05,len05)
	th05dot = thetadot!(thlen05,params)
	m05 = thlen05.mterm
	# Take the second step of RK2.
	len1 = len0 + dt*m05
	theta1 = theta0 + dt*th05dot
	# Copy the new values to a new ThetaLenType variable.
	thlen1 = new_thlen(thlen0)
	thlen1.theta = theta1
	thlen1.len = len1
	return thlen1
end
#= thetadot: Calculate the time derivative of theta;
It also calculates mterm and nterm in thlen. =#
function thetadot!(thlen::ThetaLenType, params::ParamType)
	# Extract the needed variables.
	dt, epsilon, beta = params.dt, params.epsilon, params.beta
	alpha, theta, len = thlen.alpha, thlen.theta, thlen.len
	# Calculate mterm and nterm at time 0.
	getmn!(thlen,params)
	nterm = thlen.nterm
	# Calculate the time derivative of theta.
	dth = specdiff(theta - 2*pi*alpha) + 2*pi
	d2th = specdiff(dth)
	thdot = epsilon*len^(beta-2)*d2th + nterm
	# Return thdot
	return thdot
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
function getxy(thlen::ThetaLenType)
	theta = thlen.theta
	len =  thlen.len
	xc = thlen.xc
	yc = thlen.yc
	xx,yy = getxy(theta,len,xc,yc)
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
