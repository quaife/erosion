# main.jl: The main routines to call
using Plots
using DelimitedFiles
using Statistics
using FFTW
include("basic.jl")
include("spectral.jl")
include("thetalen.jl")
include("ioroutines.jl")



# Shorthand to round a float to a number of significant digits.
sig(var::AbstractFloat, sigdig::Int=3) = round(var,sigdigits=sigdig)

# Add a variable incrementally to a jld data file.
function add_data(file::AbstractString, varlabel::AbstractString, var)
	iofile = jldopen(file, "r+")
		write(iofile, varlabel, var)
	close(iofile)
end

# Set the plot folder.
plotfolder(params::paramset) = string("../zFigs-",params.label,"/")
# If the folder exists, delete it. Then create a new folder.
function newfolder(foldername::AbstractString)
	if isdir(foldername) rm(foldername; recursive=true) end
	mkdir(foldername)
end
# readvec: Read a vector from a text file.
function readvec(filename::AbstractString)
	iostream = open(filename, "r")
	invec = readdlm(iostream, comments=true)[:,1]
	close(iostream)
	return invec
end






#--------------- MAIN ROUTINE ---------------#
# The main routine to call erosion.
function main(params::paramset)
	# Initialize the output jld file and save the parameters.
	save(params.outfile, "params", params)
	# Run the erosion simulation.
	println("BEGINNING EROSION SIMULATION")
	cputime = @elapsed(erosion(params))
	# Save the CPU time of the simulation.
	cputime = sig(cputime/3600,3)
	add_data(params.outfile, "cputime", cputime)
	println("\n\n\nCOMPLETED EROSION SIMULATION")
	println("cpu time = ", cputime, " hours.\n\n")
end

# The routine to erode a group of bodies.
function erosion(params::paramset)
	thlenden = get_thlenden(params)
	newfolder(plotfolder(params))
	nn=0; nfile = 0; tt = 0.;
	plotnsave(nfile,tt,thlenden,params)
	# Enter the time loop to apply Runge-Kutta.
	while(tt < params.tfin - 0.1*params.dt && length(thlenden.thlenvec) > 0)
		# Print statements
		nn += 1
		println("\n\n\nTIME STEP ", nn)
		println("t/tfin = ", sig(tt/params.tfin, 3))
		# Advance the variables forward one timestep with RK4.
		thlenden, dt = rungekutta2(thlenden, params)
		tt += dt
		# Plot and save the data if appropriate.
		if mod(nn, params.cntout)==0
			nfile += 1
			plotnsave(nfile,tt,thlenden,params)
		end
	end
	# Plot and save one last time with zero bodies.
	nfile += 1
	plotnsave(nfile,tt,thlenden,params)
	return thlenden,params,tt
end

# Get initial thlenden.
function get_thlenden(params::paramset)
	circdata = readvec(params.infile)
	nbods = round(Int, circdata[1])
	deleteat!(circdata,1)
	thlenvec = new_thlenvec(nbods)
	for nn = 1:nbods
		rad, xc, yc = circdata[1], circdata[2], circdata[3]
		thlenvec[nn] = circ2thlen(params.npts, rad, xc, yc)
		deleteat!(circdata,1:3)
	end
	return new_thlenden(thlenvec)
end
# Return thlen data for a circle of given radius and center.
function circ2thlen(npts::Int, rad::Float64, xc::Float64, yc::Float64)
	thlen = new_thlen()
	alpha = getalpha(npts)
	# Get the tangent angle, theta, and the total arclength, len.
	thlen.theta = 0.5*pi .+ alpha
	thlen.len = 2*pi*rad
	# Save the surface mean values, xsm and ysm.
	thlen.xsm = xc; thlen.ysm = yc
	return thlen
end





## TO DO: Other things to keep track of in addition to fixed parameters: cntout, 
###   paramdata = [label2; params.cntout; nfile; cputime]



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




