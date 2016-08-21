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
	# Create a new ThetaLenType variable.
	thlen = new_thlen()
	# alpha = s/L is the parameterization variable.
	alpha = getalpha(npts)
	# theta is the tangent angle.
	thlen.theta = 0.5*pi + 2*pi*alpha
	# len is the total arclength.
	thlen.len = 2*pi*rad
	# Save xc and yc too.
	thlen.xc = xc; thlen.yc = yc
	return thlen
end

# heaviside: The Heaviside function
function heaviside(tt)
   0.5 * (sign(tt) + 1)
end
