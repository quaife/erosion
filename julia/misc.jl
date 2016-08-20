#= misc.jl: Miscellaneous functions =#

# Create a ParamType to avoid passing all of the parameters to every function.
type ParamType
	dt::Float64; epsilon::Float64; beta::Real
end
# Function to extract the values from a ParamType variable.
function getparams(params::ParamType)
	dt = params.dt; epsilon = params.epsilon; beta = params.beta
	return dt, epsilon, beta
end

# Creat a ThetaLenType to keep all of the information of a curve.
type ThetaLenType
	npts::Integer; alpha::Vector{Float64}; theta::Vector{Float64}; len::Float64
	xc::Float64; yc::Float64; atau::Vector{Float64}; MM::Float64; NN::Vector{Float64}
end

########## Starter routines ##########
# RKstarter: Explicit second-order Runge-Kutta to start the time stepping.
function RKstarter(thlen0::ThetaLenType, params::ParamType)
	dt,epsilon,beta = getparams(params)
	theta0 = thlen0.theta
	len0 = thlen0.len


	# Get the time derivatives at t=0.
	thlen0.atau = stokes_thl_sing(theta0,len0)
	th0dot,M0,N0 = thetadot(thlen0,params)

	# Take the first half-step of RK2.
	len05 = len0 + 0.5*dt*M0
	theta05 = theta0 + 0.5*dt*th0dot


	# Save the results in a new ....
	thlen05 = deepcopy(thlen0)

	# Get the time derivatives at t=0.5*dt.
	thlen05.atau = stokes_thl_sing(theta05,len05)
	th05dot,M05,N05 = thetadot(thlen05,params)
	# Take the second step of RK2.
	len1 = len0 + dt*M05
	theta1 = theta0 + dt*th05dot

	# Remember the MM and NN values in thlen0
	thlen0.MM = M0
	thlen0.NN = N0

	# Question: will this modify thlen0 without returning it??? I think yes.

	# Copy the new values to a new ThetaLenType variable.
	thlen1 = deepcopy(thlen0)
	thlen1.theta = theta1
	thlen1.len = len1
	return thlen1
end
#= thetadot: Calculate the time derivative of theta;
Also return MM and NN while we're at it. Only used in the RKstarter. =#
function thetadot(thlen::ThetaLenType, params::ParamType)
	theta = thlen.theta
	len = thlen.len
	atau = thlen.atau

	dt,epsilon,beta = getparams(params)
	# Get the M and N terms
	MM,NN = getMN(atau,theta,len,params)
	# Calculate the time derivative of theta.
	alpha = getalpha(endof(theta))
	dth = specdiff(theta - 2*pi*alpha) + 2*pi
	d2th = specdiff(dth)
	thdot = epsilon*len^(beta-2)*d2th + NN
	return thdot, MM, NN
end

function thetadot(thlen::ThetaLenType, params::ParamType)

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
