# geos.jl
# Geometry routines

# circgeo: Creates a circle.
function circgeo(npts::Integer, rad::Float64, xsm::Float64=0.0, ysm::Float64=0.0)
	thlen = new_thlen()
	alpha = getalpha(npts)
	# Get the tangent angle, theta, and the total arclength, len.
	thlen.theta = 0.5*pi + 2*pi*alpha
	thlen.len = 2*pi*rad
	# Save the surface mean values, xsm and ysm.
	thlen.xsm = xsm; thlen.ysm = ysm
	return thlen
end
# polygongeo: Creates a polygon with number of sides, nsides.
function polygongeo(npts::Integer, nsides::Integer, 
		sigma::Float64 = 0.1, sdlen::Float64=0.5, xsm::Float64=0.0, ysm::Float64=0.0)
	thlen = new_thlen()
	alpha = getalpha(npts)
	# Get the tangent angle, theta, and smooth it.
	theta = zeros(npts)
	for nn=1:nsides
		a0 = (nn-1)/nsides
		a1 = nn/nsides
		intvl = ((a0 .<= alpha) & (alpha .<= a1))
		theta[intvl] = 0.5*pi + 2*pi*(nn-1)/nsides
	end
	theta = gaussfilter(theta - 2*pi*alpha, sigma) + 2*pi*alpha
	# Save theta, len, xsm, ysm in thlen.
	thlen.theta = theta
	thlen.len = nsides*sdlen
	thlen.xsm = xsm; thlen.ysm = ysm
	return thlen
end
# make1circ: Make a single circular geometry
function make1circ(filename::AbstractString, npts::Integer, rad::Float64, 
		xsm::Float64=0., ysm::Float64=0.)
	# Create the data vector.
	outvec = zeros(Float64,5+npts)
	outvec[1] = npts
	outvec[2] = 1		# The number of bodies.
	# Save the theta, len, xsm, ysm values in a single vector.
	alpha = getalpha(npts)
	outvec[3:2+npts] = 0.5*pi + 2*pi*alpha[:]
	outvec[3+npts] = 2*pi*rad
	outvec[4+npts] = xsm
	outvec[5+npts] = ysm
	# Write the vector to a data file.
	iostream = open(filename, "w")
	writedlm(iostream, outvec)
	close(iostream)
	return
end
# make4circs: Makes up to four circles and stores the data in a file.
function make6circs(npts::Integer, filename::AbstractString="../06circ.in")
	# Create the data vector.
	nbods = 6
	datavec = zeros(Float64, 2)
	datavec[1] = npts
	datavec[2] = nbods
	# Make some circles.
	rad = [0.39, 0.38, 0.30, 0.26, 0.25, 0.13]
	xsm,ysm = [zeros(Float64,nbods) for ii=1:2]
	xsm[1], ysm[1] = -0.40, -0.14
	xsm[2], ysm[2] = +0.38, +0.42
	xsm[3], ysm[3] = +0.60, -0.41
	xsm[4], ysm[4] = -0.22, +0.70
	xsm[5], ysm[5] = +0.07, -0.62
	xsm[6], ysm[6] = -0.73, -0.62
	# Define things.
	alpha = getalpha(npts)
	thlenvec = new_thlenvec(nbods)
	theta = 0.5*pi + 2*pi*alpha[:]
	# Save the theta, len, xsm, ysm values in a single vector.
	for nn=1:nbods
		len = 2*pi*rad[nn]
		# Create thlenvec for plotting.
		thlenvec[nn].theta = theta
		thlenvec[nn].len = len
		thlenvec[nn].xsm = xsm[nn]
		thlenvec[nn].ysm = ysm[nn]
		# Save in outvec for writing to file.
		append!(datavec, [theta; len; xsm[nn]; ysm[nn]])
	end
	writedata(datavec, filename)
	plot_curves(thlenvec, "../fig0.pdf")
	return
end
