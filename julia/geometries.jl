# geometries.jl: Construct various geometries.

# heaviside: The Heaviside function
function heaviside(tt)
   0.5 * (sign(tt) + 1)
end

# circgeo
function circgeo(npts::Integer, rad::Float64)
	# alpha = s/L is the parameterization variable.
	dalpha = 1.0/npts
	alpha = collect(range(0.5*dalpha, dalpha, npts))
	# theta is the tangent angle.
	theta = 0.5*pi - 2*pi*alpha
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
	# Parameterize the triangle by t in [0,1] and use and offset t-grid.
	dt = 1.0/npts
	tt = collect(range(0.5*dt, dt, npts))
	# alpha is the half opening angle of the triangle in radians.
	alpha = 0.5 * angle*pi/180
	# The length of the traingle's sloped face.
	len1 = 1/sin(alpha)
	# Total arclength is the sum of the 2 sloped faces and back vertical face.
	stot = 2+2*len1
	# The values of t where we change from one face to another.
	t1 = len1 / stot
	t2 = (2+len1) / stot
	# The tangent angle theta as a function of t.
	theta = 0.0*dt
	theta = alpha 
	theta += (-0.5*pi-alpha)*heaviside(tt-t1) 
	theta += (-0.5*pi-alpha)*heaviside(tt-t2)

	# Smooth the triangle.
	theta = gaussfilter(theta + 2*pi*tt, sigma) - 2*pi*tt
	# To parameterize in the CCW direction, reverse the vector and add pi.
	theta = reverse(theta) + pi
	# The point (0,xback) fixes the nose at the origin. 
	xback = 1/tan(alpha)
	return theta,stot,xback,tt
end
