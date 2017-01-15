# misc.jl

#################### Converting between x,y and theta,len ####################
# getxy!: Dispatch for input of type ThetaLenType. Only computes if they are not loaded.
function getxy!(thlen::ThetaLenType)
	# Only compute xx and yy if they are not already loaded in thlen.
	if thlen.xx==[] || thlen.yy==[]
		thlen.xx, thlen.yy = getxy(thlen.theta, thlen.len, thlen.xsm, thlen.ysm)
	end
	return
end
#= getxy: Given theta and len, reconstruct the x and y coordinates of a body.
xsm and ysm are the boundary-averaged values. =#
function getxy(theta::Vector{Float64}, len::Float64, xsm::Float64, ysm::Float64)
	# The partial derivatives dx/dalpha and dy/dalpha.
	dx = len * (cos(theta) - mean(cos(theta)))
	dy = len * (sin(theta) - mean(sin(theta)))
	# Integrate to get the x,y coordinates; result will have mean zero.
	xx = specint(dx); yy = specint(dy)
	# Move to have the correct average values.
	xx += xsm; yy += ysm
	return xx,yy
end
# getnormals: Calculate the surface normals.
function getnormals(theta::Vector{Float64})
	nx = -sin(theta)
	ny = cos(theta)
	return nx, ny 
end

#################### Geometry routines ####################
# getalpha: Calculate the parameterization variable, alpha = s/L, using an offset grid.
function getalpha(npts::Integer)
	dalpha = 1.0/npts
	return alpha = collect(range(0.5*dalpha, dalpha, npts))
end
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
function make1circ(filename::AbstractString, 
		npts::Integer, rad::Float64, xsm::Float64=0., ysm::Float64=0.)
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
	iostream = open(string(filename), "w")
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
	rad = 0.2
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
		outvec[n2-2] = 2*pi*rad
		outvec[n2-1] = xsm[nn]
		outvec[n2] = ysm[nn]
	end
	iostream = open(string(filename), "w")
	writedlm(iostream, outvec)
	close(iostream)
	return
end

#################### Plotting routines ####################
# plotcurve: Plot multiple curves from the theta-len values.
function plotcurves(thlenvec::Vector{ThetaLenType}, figname::AbstractString;
		axlims::Vector{Float64}=[3.0,1.0] )
	# Make figure of given height and preserve the aspect ratio.
	height = 400
	width = axlims[1]/axlims[2]*height
	pp = plot()
	xlim(-axlims[1],axlims[1]); ylim(-axlims[2],axlims[2])
	for ii = 1:endof(thlenvec)
		thlen = thlenvec[ii]
		if thlen.len<=0
			throw("Cannot plot a curve with non-positive length.")
		end
		getxy!(thlen)
		xx, yy = thlen.xx, thlen.yy
		pp = oplot(xx,yy,"-")
	end
	# Save the figure in a file.
	savefig(pp, figname, width=width, height=height)
	return
end
# plottheta: Plot the function theta of the first body only.
function plottheta(thlenvec::Vector{ThetaLenType}, figname::AbstractString)
	thlen = thlenvec[1]
	npts = endof(thlen.theta)
	alpha = getalpha(npts)
	pp = plot(alpha, thlen.theta, "-")
	savefig(pp, figname, width=400, height=600)
end
# plotpress: Plot the pressure distribution.
function plotpress(ytar::Vector{Float64}, ptar::Vector{Float64}, 
		figname::AbstractString, pmax::Float64=15.0)
	pp = plot(ytar, ptar, ".-")
	xlim(-1., 1.); ylim(0., pmax)
	savefig(pp, figname, width=400, height=400)
end

#################### Data IO routines ####################
# readgeodata: Reads the geometry from a data file.
# The data in the file is npts and nbods and then theta,len,xsm,yxm for each body.
function readgeodata(filename::AbstractString)
	# Open the input data file.
	iostream = open(string(filename), "r")
	invec = readdlm(iostream)
	close(iostream)
	# Extract the number of points and bodies.
	npts = round(Int,invec[1])
	nbods = round(Int,invec[2])
	# Consistency test.
	nparams = 2
	vsize = npts + 3
	if (endof(invec) != vsize*nbods+nparams)
		throw("Inconsistency in the data file."); return
	end
	# Extract the theta, len, xsm, ysm values.
	thlenvec = [new_thlen() for nn=1:nbods]
	for nn=1:nbods
		n1,n2 = n1n2(vsize,nn)
		n1 += nparams; n2 += nparams
		thlenvec[nn].theta = invec[n1:n2-3]
		testtheta(thlenvec[nn].theta)
		thlenvec[nn].len = invec[n2-2]
		thlenvec[nn].xsm = invec[n2-1]
		thlenvec[nn].ysm = invec[n2]
	end
	return thlenvec,npts
end
# testtheta: Make sure theta[end]-theta[1] = 2*pi
function testtheta(theta::Vector{Float64})
	dth = 2*pi/endof(theta)
	diff = abs(theta[end]-theta[1]-2*pi)
	if diff>2*dth
		throw("Unacceptable theta vector")
	end
end
# savexydata: Save the xy values in a data file.
function savexydata(thlenvec::Vector{ThetaLenType}, filename::AbstractString)
	nbods = endof(thlenvec)
	npts = endof(thlenvec[1].theta)
	iostream = open(string(filename), "w")
	writedlm(iostream, [npts; nbods])
	for nn=1:nbods
		getxy!(thlenvec[nn])
		xyvec = [thlenvec[nn].xx; thlenvec[nn].yy]
		writedlm(iostream, xyvec)
	end
	close(iostream)
end
