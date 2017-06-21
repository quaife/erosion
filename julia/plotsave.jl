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
	write_data(tt,thlenden,params.nouter,geomfile,densityfile)
	write_params(paramvec,paramfile)
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
function write_data(tt::Float64, thlenden::ThLenDenType, nouter::Int,
		geomfile::AbstractString, densityfile::AbstractString)
	# The quantities of interest.
	thlenvec = thlenden.thlenvec
	density = thlenden.density
	# Write the geometry data.
	iostream = open(geomfile, "w")
	npts,nbods = getnvals(thlenvec)
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
	writedlm(iostream, [tt; npts; nbods; nouter; density])
	close(iostream)
end
# write_params: Write the important parameters in an output file.
function write_params(paramvec::Array, filename::AbstractString)
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
# readvec: Read a vector from a data file.
function readvec(filename::AbstractString)
	iostream = open(filename, "r")
	invec = readdlm(iostream)[:,1]
	close(iostream)
	return invec
end
# readthlenfile: Reads the geometry from a data file.
# The data in the file is npts and nbods and then theta,len,xsm,yxm for each body.
function read_thlen_file(filename::AbstractString)
	# Read the data file.
	invec = readvec(filename)
	# Extract the number of points and bodies.
	npts = round(Int,invec[1])
	nbods = round(Int,invec[2])
	deleteat!(invec,1:2)
	# Consistency test.
	assert( endof(invec) == nbods*(npts+3))
	# Extract thlenvec for each body.
	thlenvec = new_thlenvec(nbods)
	for nn=1:nbods
		thlenvec[nn].theta = invec[1:npts]
		thlenvec[nn].len = invec[npts+1]
		thlenvec[nn].xsm = invec[npts+2]
		thlenvec[nn].ysm = invec[npts+3]
		test_theta(thlenvec[nn].theta)
		deleteat!(invec,1:npts+3)
	end
	return thlenvec
end
# geom2thlen: Convert a geom.dat file to thlen.in file.
function geom2thlen(infile::AbstractString, outfile::AbstractString)
	# Read the data file.
	invec = readvec(infile)
	# Extract the parameters.
	tt = invec[1]
	npts = round(Int,invec[2])
	nbods = round(Int,invec[3])
	deleteat!(invec,1:3)
	# Consistency test.
	assert( endof(invec) == nbods*(3*npts+3))
	# Start to write data.
	iostream = open(outfile, "w")
	writedlm(iostream, [npts; nbods])
	# Write theta, len, xsm, ysm to outfile and skip over xx and yy.
	for nn=1:nbods
		theta = invec[1:npts]
		len = invec[npts+1]
		xsm = invec[npts+2]
		ysm = invec[npts+3]
		test_theta(theta)
		writedlm(iostream, [theta; len; xsm; ysm])
		# Delete theta, len, xsm, ysm
		deleteat!(invec,1:npts+3)
		# Delete xx and yy too.
		deleteat!(invec,1:2*npts)
	end
	close(iostream)
end
