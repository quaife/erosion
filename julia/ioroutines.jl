# misc.jl
# IO routines for plotting and saving data.

# plotnsave: Calls plotcurves() and savedata()
function plotnsave(nfile::Int, tt::Float64, thlenden::ThLenDenType, params::ParamType)
	# Preliminary stuff.
	println("\n\n\nOUTPUT NUMBER ", nfile)
	# The file names.
	datafolder, plotfolder = getfoldernames()
	nfilestr = lpad(nfile,4,0)
	geomfile = string(datafolder,"geom",nfilestr,".dat")
	densityfile = string(datafolder,"density",nfilestr,".dat")
	paramsfile = string(datafolder,"params.dat")
	# Plot the shapes.
	plotfile = string(plotfolder,"shape",nfilestr,".pdf")
	plot_curves(thlenden.thlenvec,plotfile)
	# Compute the density functions.
	getstress!(thlenden, params)
	compute_denrot!(thlenden, params)
	# Write the data to a file.
	save_geo_density(tt,thlenden,geomfile,densityfile)
	save_params(params,nfile,paramsfile)
	return
end

#--------------- WRITING DATA ---------------#
# writedata: Write generic data to a file.
function writedata(data::Vector, filename::AbstractString)
	iostream = open(filename, "w")
	writedlm(iostream, data)
	close(iostream)
	return
end
#= save_geo_density: Save the geometry data (theta,len,xsm,ysm,xx,yy) 
and the density-function data in a file. =#
function save_geo_density(tt::Float64, thlenden::ThLenDenType,
		geomfile::AbstractString, densityfile::AbstractString)
	# Write the geometry data.
	thlenvec = thlenden.thlenvec
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
	densitydata = [thlenden.density; thlenden.denrot]
	writedata(densitydata,densityfile)
	return
end
# save_params: Write the important parameters in an output file.
function save_params(params::ParamType, nfile::Int, outparamsfile::AbstractString)
	paramvec = readvec("params.dat")
	cputime = round( (time()-params.cput0)/60. , 2)
	label1 = "# Input Parameters: geoinfile, nouter, tfin, dtout, dtfac, epsfac, sigfac, iffm, fixarea"
	label2 = "# Calculated Parameters: npts, cntout, last file number, cputime (minutes)"
	paramdata = [label1; paramvec; label2; 
		params.npts; params.cntout; nfile; cputime]
	writedata(paramdata, outparamsfile)
	return
end
#--------------- HANDLING FOLDERS ---------------#
# getfoldernames: Set the name of the data and plot folders.
function getfoldernames()
	datafolder = "../datafiles/run/"
	plotfolder = "../figs/"
	return datafolder, plotfolder
end
# newfolder: If the folder exists, delete it. Then create a new folder.
function newfolder(foldername::AbstractString)
	if isdir(foldername)
		rm(foldername; recursive=true)
	end
	mkdir(foldername)
	return
end

#--------------- READING DATA ---------------#
# read_thlen_file: Read a thlen file.
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
		test_theta_ends(thlenvec[nn].theta)
		test_theta_means(thlenvec[nn].theta)
		deleteat!(invec,1:npts+3)
	end
	return thlenvec
end
# read_geom_file: Read a geometry file.
#= The data in the file is t (physical time), npts, nbods, 
then theta, len, xsm, ysm, xx, yy for each body. =#
function read_geom_file(filename::AbstractString)
	# Read the data file.
	invec = readvec(filename)
	# Extract the parameters.
	tt = invec[1]
	npts = round(Int,invec[2])
	nbods = round(Int,invec[3])
	deleteat!(invec,1:3)
	# Consistency test.
	assert( endof(invec) == nbods*(3*npts+3))
	# Read theta, len, xsm, ysm, xx, yy and save in thlenvec.
	thlenvec = new_thlenvec(nbods)
	for nn=1:nbods
		# Read theta, len, xsm, ysm.
		thlenvec[nn].theta = invec[1:npts]
		thlenvec[nn].len = invec[npts+1]
		thlenvec[nn].xsm = invec[npts+2]
		thlenvec[nn].ysm = invec[npts+3]
		test_theta_ends(thlenvec[nn].theta)
		test_theta_means(thlenvec[nn].theta)
		# Delete theta, len, xsm, ysm
		deleteat!(invec,1:npts+3)
		# Read xx and yy.
		thlenvec[nn].xx = invec[1:npts]
		thlenvec[nn].yy = invec[npts+1:2*npts]
		# Delete xx and yy.
		deleteat!(invec,1:2*npts)
	end
	return tt,thlenvec
end
# read_density_file: Read the density file to get the density and rotated density.
function read_density_file(filename::AbstractString)
	dendata = readvec(filename)
	nn = endof(dendata)
	assert(iseven(nn))
	nhalf = div(nn,2)
	density = dendata[1:nhalf]
	denrot = dendata[nhalf+1:nn]
	return density, denrot
end
# readvec: Read a vector from a data file.
function readvec(filename::AbstractString)
	iostream = open(filename, "r")
	invec = readdlm(iostream)[:,1]
	close(iostream)
	return invec
end

#--------------- OTHER ---------------#
# geom2thlen: Convert a geom.dat file to thlen.in file. =#
function geom2thlen(geofile::AbstractString, thlenfile::AbstractString)
	tt,thlenvec = read_geom_file(geofile)
	npts,nbods = getnvals(thlenvec)
	datavec = zeros(Float64, 2)
	datavec[1] = npts
	datavec[2] = nbods
	for nn=1:nbods
		test_theta_ends(thlenvec[nn].theta)
		test_theta_means(thlenvec[nn].theta)
		append!(datavec, [thlenvec[nn].theta; 
			thlenvec[nn].len; thlenvec[nn].xsm; thlenvec[nn].ysm])
	end
	writedata(datavec, thlenfile)
end
# Generic version that automatically names the files.
function geom2thlen(geomnum::AbstractString)
	geofile = string("../datafiles/run/geom", geomnum, ".dat")
	thlenfile = "../thlen.in"
	geom2thlen(geofile,thlenfile)
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
