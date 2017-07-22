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
	save_geo_density(tt,thlenden,geomfile,densityfile)
	save_params(paramvec,cnt,paramfile)
	# Plot the shapes.
	plotfile = string(plotfolder,"shape",cntstr,".pdf")
	plot_curves(thlenden.thlenvec,plotfile)
	# Plot the pressures.
	#pressfile = string(plotfolder,"pressure",cntstr,".pdf")
	#pressure = compute_pressure(thlenden,params.nouter)
	#plot_pressure(pressure,pressfile)
	return
end

#--------------- SAVING DATA ---------------#
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
# write_params: Write the important parameters in an output file.
function save_params(paramvec::Array, cnt::Int, filename::AbstractString)
	label1 = "# Input Parameters: geoinfile, nouter, tfin, dtout, dtfac, epsfac, sigfac, iffm, fixarea"
	label2 = "# Calculated Parameters: dtoutexact, cntout, cputime (minutes), lastcnt"
	paramdata = [label1; paramvec[1:end-3]; 
		label2; paramvec[end-2:end-1]; round(paramvec[end],2); cnt]
	writedata(paramdata, filename)
	return
end
# writedata: Write generic data to a file.
function writedata(data::Vector, filename::AbstractString)
	iostream = open(filename, "w")
	writedlm(iostream, data)
	close(iostream)
	return
end
# newfolder: If the folder exists, delete it and create a new one.
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
		test_theta(thlenvec[nn].theta)
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
		test_theta(thlenvec[nn].theta)
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
	nden = endof(dendata)
	assert(iseven(nden))
	density = dendata[1:nden/2]
	denrot = dendata[nden/2+1:nden]
	return density, denrot
end
# readvec: Read a vector from a data file.
function readvec(filename::AbstractString)
	iostream = open(filename, "r")
	invec = readdlm(iostream)[:,1]
	close(iostream)
	return invec
end
#= REDO using read_geom_file
# geom2thlen: Convert a geom.dat file to thlen.in file. =#


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
	return
end
