# Test the curvature driven flow

# Curvature-driven flow.
function cdf(npts::Int, dt::Float64, epsilon::Float64,
		tfin::Float64=0.1, len0::Float64=2*pi*0.2, lenevo::Int=0)
	# Calculate parameters and initial shape.
	nsteps = round(Int,tfin/dt)
	params = ParamType(dt,epsilon,0.,0.,lenevo,0)
	thlenvec0 = makeshape(npts,len0,a1)
	# Start the simulation with RKstarter.
	thlenvec1 = RKstarter!(thlenvec0,params,noatau)

	# Continue the simulation with the multi-step method.
	for nn = 2:nsteps
		advance_thetalen!(thlenvec1,thlenvec0,params)
	end
	return thlenvec1
end

# Make the initial shape.
function makeshape(npts::Int, len::Float64, a1::Float64, 
		xsm::Float64=0., ysm::Float64=0.)
	alpha = getalpha(npts)
	# theta for a circle.
	theta = 0.5*pi + 2*pi*alpha[:]
	# Add a sin component.
	theta += a1*sin(2*pi*alpha[:])
	# Put the shape in a thlenvec.
	thlen = new_thlen()
	thlen.theta = theta
	thlen.len = len
	thlen.xsm = xsm
	thlen.ysm = ysm
	thlenvec = [thlen]
	return thlenvec
end

# Function to pass into RKstarter without computing abs tau.
function noatau(thlenv::Vector{ThetaLenType}, params::ParamType)
	return
end

