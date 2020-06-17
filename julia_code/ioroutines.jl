# IO routines for plotting and saving data.







# plotnsave
function plotnsave(nfile::Int, tt::Float64, thlenden::ThLenDenType, params::ParamType)
	# Preliminary stuff.
	println("\n\n\nOUTPUT NUMBER ", nfile)
	# The file names.
	plotfolder = string("../zFigs", run_label)


	# Plot the shapes.
	plotfile = string(plotfolder,"shape",nfilestr,".pdf")
	plot_curves(thlenden.thlenvec,plotfile)
	# Compute the density functions.
	getstress!(thlenden, params)
	compute_density!(thlenden, params, rotation=true)



	# Write the data to a file.

	iostream = jldopen(datafile, "r+")
		write(iostream, "y", y)
	close(iostream)

	## save_geo_density(tt,thlenden,geomfile,densityfile)
	## save_pinfo(params,nfile,pinfofile)

	return
end

#params.paramsfile

# getfoldernames: Set the name of the data and plot folders.
function getfoldernames(paramsfile::AbstractString)
	datafolder = string("../output_data/run_",paramsfile,"/")
	plotfolder = string("../zFigsErosion_",paramsfile,"/")
	return datafolder, plotfolder
end

###   paramdata = [label2; params.cntout; nfile; cputime]







#--------------- HANDLING FOLDERS ---------------#
# newfolder: If the folder exists, delete it. Then create a new folder.
function newfolder(foldername::AbstractString)
	if isdir(foldername) rm(foldername; recursive=true) end
	mkdir(foldername)
end





#--------------- PLOTTING DATA ---------------#
# plotcurve: Plot multiple curves from the theta-len values.
function plot_curves(thlenvec::Vector{ThetaLenType}, figname::AbstractString)	
	# Make figure of given height and preserve the aspect ratio.
	axlims = [1.0,1.0]
	height = 400
	width = axlims[1]/axlims[2]*height
	plt = plot(xlim=(-axlims[1],axlims[1]), ylim=(-axlims[2],axlims[2]), size=(width,height),leg=false)
	for ii = 1:lastindex(thlenvec)
		thlen = thlenvec[ii]
		if thlen.len<=0
			throw("Cannot plot a curve with non-positive length.")
		end
		getxy!(thlen)
		xx, yy = thlen.xx, thlen.yy
		plot!(plt, xx,yy,color="black")
	end
	# Save the figure in a file.
	savefig(plt, figname)
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
	@assert length(invec) == nbods*(npts+3)
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

# readvec: Read a vector from a data file.
function readvec(filename::AbstractString)
	iostream = open(filename, "r")
	invec = readdlm(iostream, comments=true)[:,1]
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
	geofile = string("../output_data/run/geom", geomnum, ".dat")
	thlenfile = "../thlen.in"
	geom2thlen(geofile,thlenfile)
end



