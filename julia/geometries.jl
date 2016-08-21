# geometries.jl: Construct various geometries.

# getalpha: Calculate the parameterization variable, alpha = s/L.
function getalpha(npts::Integer)
    # alpha = s/L is the parameterization variable.
    dalpha = 1.0/npts
    # Use an offset grid.
    alpha = collect(range(0.5*dalpha, dalpha, npts))
    return alpha
end

# circgeo: Creates a circle.
function circgeo(npts::Integer, rad::Float64, xc::Float64=0.0, yc::Float64=0.0)
    # alpha = s/L is the parameterization variable.
    alpha = getalpha(npts)
	# theta is the tangent angle.
	theta = 0.5*pi + 2*pi*alpha
	# len is the total arclength.
	len = 2*pi*rad
	# Save the data in a ThetaLenType variable.
    thlen = new_thlen()
	thlen.theta = theta
    thlen.len = len
    thlen.xc = xc
    thlen.yc = yc
	return thlen
end

#= trigeo: Construct a triangle geometry within the theta-L framework.
npts: The number of points on the triangle.
angle: The triangle's front opening angle in degrees.
sigma: The amount of smoothing. 
The triangle has the following properties:
- It is up-down symmetric.
- It's nose is at the origin.
- It has a vetical extent of 2 units.
- It is parameterized by t, i.e. (x(t), y(t)), for t in [0,1]
- It is parameterized in the CCW direction. =#
function trigeo(npts::Integer, angle, sigma)
	#= Initially, the triangle is parameterized in CW direction, 
	but that will be reversed at the end. =#
	# Parameterize the triangle by tt in [0,1] and use and offset tt-grid.
	tt = getalpha(npts)
	# hang is the half opening angle of the triangle in radians.
	hang = 0.5*angle*pi/180
	# The length of the traingle's sloped face.
	len1 = 1/sin(hang)
	# Total arclength is the sum of the 2 sloped faces and back vertical face.
	stot = 2+2*len1
	# The values of t where we change from one face to another.
	t1 = len1 / stot
	t2 = (2+len1) / stot
	# The tangent angle theta as a function of t.
	theta = 0.0*dt
	theta = hang 
	theta += (-0.5*pi-hang)*heaviside(tt-t1) 
	theta += (-0.5*pi-hang)*heaviside(tt-t2)

	# Smooth the triangle.
	theta = gaussfilter(theta + 2*pi*tt, sigma) - 2*pi*tt
	# To parameterize in the CCW direction, reverse the vector and add pi.
	theta = reverse(theta) + pi
	# The point (0,xback) fixes the nose at the origin. 
	xback = 1/tan(hang)
	return theta, stot, xback
end

# heaviside: The Heaviside function
function heaviside(tt)
   0.5 * (sign(tt) + 1)
end
