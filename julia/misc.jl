# misc.jl

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
	iostream = open(filename, "w")
	writedlm(iostream, outvec)
	close(iostream)
	return
end

#################### Plotting routines ####################
# plotcurve: Plot multiple curves from the theta-len values.
function plotcurves(thlenvec::Vector{ThetaLenType}, figname::AbstractString)	
	# Make figure of given height and preserve the aspect ratio.
	axlims = [1.0,1.0]
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

#################### Data IO routines ####################
# readthlenfile: Reads the geometry from a data file.
# The data in the file is npts and nbods and then theta,len,xsm,yxm for each body.
function readthlenfile(filename::AbstractString)
	# Open the input data file.
	iostream = open(filename, "r")
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
	return thlenvec
end
# testtheta: Test that the theta vector is reasonable.
function testtheta(theta::Vector{Float64})
	npts = endof(theta)
	# 1) Make sure that the jump between the endpoints is 2*pi.
	# Linear extrapolation to estimate theta at alpha=0 from both sides.
	th0left = 1.5*theta[end] - 0.5*theta[end-1] - 2*pi
	th0right = 1.5*theta[1] - 0.5*theta[2]
	# Compare the two extrpaolations.
	th0diff = abs(th0left - th0right)
	thresh = 0.2
	if th0diff > thresh
		throw(string("Unacceptable theta vector, ", 
			"the endpoints do not match: ", signif(th0diff,3), " > ", signif(thresh,3) ))
	end
	# 2) Make sure that cos(theta) and sin(theta) have zero mean.
	m1 = mean(cos(theta))
	m2 = mean(sin(theta))
	maxmean = maximum(abs([m1,m2]))
	thresh = 10./npts
	if maxmean > thresh
		throw(string("Unacceptable theta vector, ",
			"the means are not right: ", signif(maxmean,3), " > ", signif(thresh,3) ))
	end
end
# savedata: Save the all of the data (theta,len,xsm,ysm,xx,yy) in a file.
function savedata(thlenvec::Vector{ThetaLenType}, tt::Float64, filename::AbstractString)
	nbods = endof(thlenvec)
	npts = endof(thlenvec[1].theta)
	iostream = open(filename, "w")
	writedlm(iostream, [tt; npts; nbods])
	for nn=1:nbods
		getxy!(thlenvec[nn])
		datavec = [thlenvec[nn].theta; thlenvec[nn].len; 
			thlenvec[nn].xsm; thlenvec[nn].ysm; 
			thlenvec[nn].xx; thlenvec[nn].yy]
		writedlm(iostream, datavec)
	end
	close(iostream)
end
# plotnsave: Calls plotcurves() and savedata()
function plotnsave(thlenvec::Vector{ThetaLenType}, datafolder::AbstractString, 
		plotfolder::AbstractString, tt::Float64, cnt::Integer)
	# Save the data.
	cntstr = lpad(cnt,4,0)
	savefile = string(datafolder,"geom",cntstr,".dat")
	savedata(thlenvec,tt,savefile)
	# Plot the shapes.
	plotfile = string(plotfolder,"shape",cntstr,".pdf")
	plotcurves(thlenvec,plotfile)
end
# newfolder: If the folder exists, delete it and create a new one.
function newfolder(foldername::AbstractString)
	if isdir(foldername)
		rm(foldername; recursive=true)
	end
	mkdir(foldername)
end
# paramsout: Save important parameters in an output file.
function writeparams(filename::AbstractString, paramvec::Array)
	label1 = "# Input Parameters: geoinfile, nouter, tfin, dtout, dtfac, epsfac, sigfac, iffm, fixarea"
	label2 = "# Calculated Parameters: dtoutexact, cntout, cputime (minutes)"
	writevec = [label1; paramvec[1:end-3]; label2; paramvec[end-2:end-1]; round(paramvec[end],2)]
	iostream = open(filename, "w")
	writedlm(iostream, writevec)
	close(iostream)
end
