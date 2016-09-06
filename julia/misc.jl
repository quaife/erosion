# misc.jl

#################### Geometry routines ####################
# getalpha: Calculate the parameterization variable, alpha = s/L.
function getalpha(npts::Integer)
	# alpha = s/L is the parameterization variable.
	dalpha = 1.0/npts
	# Use an offset grid.
	alpha = collect(range(0.5*dalpha, dalpha, npts))
	return alpha
end
# circgeo: Creates a circle.
function circgeo(npts::Integer, rad::Float64, xsm::Float64=0.0, ysm::Float64=0.0)
	# Create a new ThetaLenType variable.
	thlen = new_thlen()
	# alpha = s/L is the parameterization variable.
	alpha = getalpha(npts)
	# theta is the tangent angle.
	thlen.theta = 0.5*pi + 2*pi*alpha
	# len is the total arclength.
	thlen.len = 2*pi*rad
	# Save xsm and ysm too.
	thlen.xsm = xsm; thlen.ysm = ysm
	return thlen
end
# polygongeo: Creates a polygon with number of sides, nsides.
function polygongeo(npts::Integer, nsides::Integer, 
		sigma::Float64 = 0.1, sdlen::Float64=0.5, xsm::Float64=0.0, ysm::Float64=0.0)
	# Create a new ThetaLenType variable.
	thlen = new_thlen()
	# alpha = s/L is the parameterization variable.
	alpha = getalpha(npts)
	# theta is the tangent angle.
	theta = zeros(npts)
	for nn=1:nsides
		a0 = (nn-1)/nsides
		a1 = nn/nsides
		intvl = ((a0 .<= alpha) & (alpha .<= a1))
		theta[intvl] = 0.5*pi + 2*pi*(nn-1)/nsides
	end
	# Smooth theta
	theta = gaussfilter(theta - 2*pi*alpha, sigma) + 2*pi*alpha
	# Save theta in the ThetaLenType object.
	thlen.theta = theta
	# len is the total arclength.
	thlen.len = nsides*sdlen
	# Save xsm and ysm too.
	thlen.xsm = xsm; thlen.ysm = ysm
	return thlen
end
##################################################


#################### Other routines ####################
#= getxy: Given theta and len, reconstruct the x and y coordinates of a body.
xsm and ysm are the boundary-averaged values.
While we're at it, we can also calculate the normal direcations. =#
function getxy(theta::Vector{Float64}, len::Float64, xsm::Float64, ysm::Float64)
	# The increments of dx and dy
	dx = len * (cos(theta) - mean(cos(theta)))
	dy = len * (sin(theta) - mean(sin(theta)))
	# Integrate to get the x,y coordinates; result will have mean zero.
	xx = specint(dx)
	yy = specint(dy)
	# Move to have the correct average values.
	xx += xsm
	yy += ysm
	return xx,yy
end
# getxy!: Dispatch for input of type ThetaLenType. Only computes if they are not loaded.
function getxy!(thlen::ThetaLenType)
	# Only compute xx and yy if they are not already loaded in thlen.
	if thlen.xx==[] || thlen.yy==[]
		thlen.xx, thlen.yy = getxy(thlen.theta, thlen.len, thlen.xsm, thlen.ysm)
	end
	return
end
# getnormals: Calculate the surface normals.
function getnormals(theta::Vector{Float64})
	nx = -sin(theta)
	ny = cos(theta)
	return nx, ny 
end
# plotcurve: Plot a curve from the theta-len values along with the initial shape.
function plotsinglecurve!(thlen1::ThetaLenType, thlen0::ThetaLenType, cnt::Integer; axlim::Real=0.5)
	# Compute the xy coordinates if they are not already loaded in thlen.
	getxy!(thlen0); getxy!(thlen1)
	x0, y0 = thlen0.xx, thlen0.yy
	x1, y1 = thlen1.xx, thlen1.yy
	# Plot the curves.
	p1 = plot(x0,y0,"-", x1,y1,"-")
	xlim(-axlim,axlim); ylim(-axlim,axlim)
	# Save the figures in a folder.
	figname = string("../figs/fig",string(cnt),".pdf")
	savefig(p1, figname, width=500, height=500)
	return
end
# plotcurve: Plot multiple curves from the theta-len values.
function plotcurves!(thlenvec::Vector{ThetaLenType}, cnt::Integer; axlim::Real=0.5)
	p1 = plot()
	xlim(-axlim,axlim); ylim(-axlim,axlim)
	for ii = 1:endof(thlenvec)
		thlen = thlenvec[ii]
		# Compute the xy coordinates if they are not already loaded in thlen.
		getxy!(thlen)
		xx, yy = thlen.xx, thlen.yy
		# Plot the curves.
		p1 = oplot(xx,yy,"-")
	end
	# Save the figures in a folder.
	figname = string("../figs/fig",string(cnt),".pdf")
	savefig(p1, figname, width=500, height=500)
	return
end
############################################################


#################### Starter routines ####################
#= festep: Take a single forward Euler step of theta, len, xsm, and ysm.
Note: Do not use the dt inside params because I might want to input 
something else, like 0.5*dt for the Runge-Kutta starter. =#
function festep(dt::Float64, thdot::Vector{Float64}, 
		theta0::Vector{Float64}, len0::Float64, xsm0::Float64, ysm0::Float64,
		ldot::Float64, xsmdot::Float64, ysmdot::Float64)
	theta1 = theta0 + dt*thdot
	len1 = len0 + dt*ldot
	xsm1 = xsm0 + dt*xsmdot
	ysm1 = ysm0 + dt*ysmdot
	return theta1, len1, xsm1, ysm1
end
# festep: Dispatch for ThetaLenType.
# Allow the starting point and the point where the derivatives are taken to be different.
function festep(dt::Float64, thdot::Vector{Float64}, 
		thlen0::ThetaLenType, thlendots::ThetaLenType)
	thlen1 = new_thlen()
	thlen1.theta, thlen1.len, thlen1.xsm, thlen1.ysm = festep(dt, thdot, 
		thlen0.theta, thlen0.len, thlen0.xsm, thlen0.ysm, 
		thlendots.mterm, thlendots.xsmdot, thlendots.ysmdot)
	return thlen1
end

#= RKstarter!: Explicit second-order Runge-Kutta to start the time stepping.
It also calculates mterm, nterm, xsmdot, ysmdot and saves them in thlen0. =#
function RKstarter!(thlen0::ThetaLenType, params::ParamType)
	dt = params.dt
	# Compute the stress at t=0.
	stokes!([thlen0])
	# Calculate the time derivatives: thdot, mterm, xsmdot, ysmdot.
	th0dot = thetadot!(thlen0,params)
	# Take the first half-step of RK2.
	thlen05 = festep(0.5*dt, th0dot, thlen0, thlen0)
	# Compute the stress at t=0.5*dt.
	stokes!([thlen05])
	# Calculate the time derivatives: thdot, mterm, xsmdot, ysmdot.
	th05dot= thetadot!(thlen05,params)
	# Take the second half-step of RK2 where the starting point is thlen0 but
	# the derivatives are calculated at thlen05.
	thlen1 = festep(dt, th0dot, thlen0, thlen05)
	# Return.
	return thlen1
end
# RKstarter!: Dispatch for vector of ThetaLenType to handle multiple bodies.
function RKstarter!(thlenvec0::Vector{ThetaLenType}, params::ParamType)
	dt = params.dt
	nbods = endof(thlenvec0) 
	# Initialize vectors of ThetaLenType.
	thlenvec05 = Array(ThetaLenType, nbods)
	thlenvec1 = Array(ThetaLenType, nbods)
	# Compute the stress at t=0.
	stokes!(thlenvec0)	
	# For each body, take the first step of RK2.
	for ii = 1:nbods
		# Need thlen0 for each body.
		thlen0 = thlenvec0[ii]
		# Calculate the time derivatives: thdot, mterm, xsmdot, ysmdot.
		thdot = thetadot!(thlen0,params)
		# Take the first step of RK2.
		thlenvec05[ii] = festep(0.5*dt, thdot, thlen0, thlen0)		
	end
	# Compute the stress at t=0.5*dt.
	stokes!(thlenvec05)	
	# For each body, take the second step of RK2.
	for ii = 1:nbods
		# Need both thlen0 and thlen05 for each body.
		thlen0 = thlenvec0[ii]
		thlen05 = thlenvec05[ii]
		# Calculate the time derivatives: thdot, mterm, xsmdot, ysmdot.
		thdot = thetadot!(thlen05,params)
		# Take the second step of RK2.
		thlenvec1[ii] = festep(dt, thdot, thlen0, thlen05)		
	end
	return thlenvec1
end
############################################################