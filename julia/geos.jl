# geos.jl

#################### Geometry routines ####################
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
function make4circs(filename::AbstractString, npts::Integer, nbods::Integer)
	# Create the data vector.
	nparams = 2
	vsize = npts + 3
	outvec = zeros(Float64,vsize*nbods+nparams)
	outvec[1] = npts
	outvec[2] = nbods
	# For now, make some circles.
	rad = [0.2, 0.2, 0.02, 0.1]
	xsm,ysm = [zeros(Float64,4) for ii=1:2]
	xsm[1], ysm[1] = +0.1, +0.4
	xsm[2], ysm[2] = -0.0, -0.4
	xsm[3], ysm[3] = +0.4, +0.1
	xsm[4], ysm[4] = -0.1, +0.0
	alpha = getalpha(npts)
	# Save the theta, len, xsm, ysm values in a single vector.
	for nn=1:nbods
		n1 = vsize*(nn-1)+1 + nparams
		n2 = vsize*nn + nparams
		outvec[n1:n2-3] = 0.5*pi + 2*pi*alpha[:]
		outvec[n2-2] = 2*pi*rad[nn]
		outvec[n2-1] = xsm[nn]
		outvec[n2] = ysm[nn]
	end
	iostream = open(filename, "w")
	writedlm(iostream, outvec)
	close(iostream)
	return
end
