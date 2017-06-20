# misc.jl
# IO routines for plotting and saving data.

# plotnsave: Calls plotcurves() and savedata()
function plotnsave(thlenden::ThLenDenType, params::ParamType, paramvec::Vector,
		datafolder::AbstractString, plotfolder::AbstractString,
		tt::Float64, cnt::Integer)
	# The file names.
	cntstr = lpad(cnt,4,0)
	geomfile = string(datafolder,"geom",cntstr,".dat")
	densityfile = 	string(datafolder,"density",cntstr,".dat")
	paramfile = string(datafolder,"params.dat")
	# Write the data to a file.
	write_data(tt,thlenden,geomfile,densityfile)
	write_params(paramfile,paramvec)
	# Plot the shapes.
	plotfile = string(plotfolder,"shape",cntstr,".pdf")
	pressfile = string(plotfolder,"pressure",cntstr,".pdf")
	pressure = getpressure(thlenden,params)
	plot_curves(thlenden.thlenvec,plotfile)
	#plot_pressure(pressure,pressfile)
end

#--------------- SAVING DATA ---------------#
#= write_data: Write the geometry data (theta,len,xsm,ysm,xx,yy) 
and the density-function data in a file. =#
function write_data(tt::Float64, thlenden::ThLenDenType, 
		geomfile::AbstractString, densityfile::AbstractString)
	# The quantities of interest.
	thlenvec = thlenden.thlenvec
	density = thlenden.density
	# Write the geometry data.
	iostream = open(geomfile, "w")
	nbods = endof(thlenvec)
	npts = endof(thlenvec[1].theta)
	label = "# Parameters (time, npts, nbods), then geometry (theta, len, xsm, ysm, x, y) for each body"
	writedlm(iostream, [label; tt; npts; nbods])
	for nn=1:nbods
		getxy!(thlenvec[nn])
		datavec = [thlenvec[nn].theta; thlenvec[nn].len; 
			thlenvec[nn].xsm; thlenvec[nn].ysm; 
			thlenvec[nn].xx; thlenvec[nn].yy]
		writedlm(iostream, datavec)
	end
	close(iostream)
	# Write the density data.
	iostream = open(densityfile, "w")
	writedlm(iostream, [tt; density])
	close(iostream)
end
# write_params: Write the important parameters in an output file.
function write_params(filename::AbstractString, paramvec::Array)
	label1 = "# Input Parameters: geoinfile, nouter, tfin, dtout, dtfac, epsfac, sigfac, iffm, fixarea"
	label2 = "# Calculated Parameters: dtoutexact, cntout, cputime (minutes)"
	writevec = [label1; paramvec[1:end-3]; label2; paramvec[end-2:end-1]; round(paramvec[end],2)]
	iostream = open(filename, "w")
	writedlm(iostream, writevec)
	close(iostream)
end
# newfolder: If the folder exists, delete it and create a new one.
function newfolder(foldername::AbstractString)
	if isdir(foldername)
		rm(foldername; recursive=true)
	end
	mkdir(foldername)
end

#--------------- PLOTTING DATA ---------------#
# plotcurve: Plot multiple curves from the theta-len values.
function plot_curves(thlenvec::Vector{ThetaLenType}, figname::AbstractString)	
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
function plot_pressure(pressure::Vector{Float64}, figname::AbstractString)
	# Make figure of given height and preserve the aspect ratio.
	height = 400
	width = 600
	pp = plot()
	# Plot the pressure
	pp = plot(pressure)
	savefig(pp, figname, width=width, height=height)
end

#--------------- READING DATA ---------------#
# readthlenfile: Reads the geometry from a data file.
# The data in the file is npts and nbods and then theta,len,xsm,yxm for each body.
function read_thlen_file(filename::AbstractString)
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
		test_theta(thlenvec[nn].theta)
		thlenvec[nn].len = invec[n2-2]
		thlenvec[nn].xsm = invec[n2-1]
		thlenvec[nn].ysm = invec[n2]
	end
	return thlenvec
end
# testtheta: Test that the theta vector is reasonable.
function test_theta(theta::Vector{Float64})
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
	thresh = 20./npts
	if maxmean > thresh
		throw(string("Unacceptable theta vector, ",
			"the means are not right: ", signif(maxmean,3), " > ", signif(thresh,3) ))
	end
end
