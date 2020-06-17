# main.jl: The main routines to call
using JLD
using DelimitedFiles
using Plots
using Statistics

include("basic.jl")
include("spectral.jl")
include("thetalen.jl")
 
# Add a variable incrementally to a jld data file.
function add_data(filename::AbstractString, varlabel::AbstractString, var)
	jldopen(filename, "r+") do file
	file[varlabel] = var
	#write(file,varlabel,var)
	end
end

#--------------- MAIN ROUTINE ---------------#
# The main routine to call erosion.
function main(params::ParamSet)
	# Initialize the output jld file and save the parameters.
	save(params.outfile, "params", params)
	# Run the erosion simulation.
	println("\nBEGINNING EROSION SIMULATION")
	cputime = @elapsed(erosion(params))
	# Save the CPU time of the simulation.
	cpu_hours = round(cputime/3600, sigdigits=3)
	add_data(params.outfile, "cpu_hours", cpu_hours)
	println("\n\n\nCOMPLETED EROSION SIMULATION")
	println("cpu time = ", cpu_hours, " hours.\n\n")
end

# The routine to erode a group of bodies.
function erosion(params::ParamSet)
	# Initialize.
	thlenden = get_thlenden(params)
	newfolder(plotfolder(params))
	nn=0; nout = 0; tt = 0.0;
	# Enter the time loop to apply Runge-Kutta.
	while(tt < params.tfin && length(thlenden.thlenvec) > 0)
		# Plot and save the data if appropriate.
		if mod(nn, params.outstride) == 0
			plotnsave(thlenden,params,nout,tt)
			nout += 1
		end
		# Advance the variables forward one timestep with RK4.
		nn += 1
		println("\n\nTIME STEP ", nn)
		thlenden, dt = rungekutta2(thlenden, params)
		tt += dt
	end
	# Plot and save one last time with zero bodies.
	plotnsave(thlenden,params,nout,tt)
	add_data(params.outfile, "noutputs", nout)
end

# Initialize thlenden from the input circle file.
function get_thlenden(params::ParamSet)
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
# readvec: Read a vector from a text file.
function readvec(filename::AbstractString)
	iostream = open(filename, "r")
	invec = readdlm(iostream, comments=true)[:,1]
	close(iostream)
	return invec
end
# Convert the circle data to thlen data.
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

# Set the plot folder.
plotfolder(params::ParamSet) = string("../zFigs-",params.label,"/")
# If the folder exists, delete it. Then create a new folder.
function newfolder(foldername::AbstractString)
	if isdir(foldername) rm(foldername; recursive=true) end
	mkdir(foldername)
end

# Plot and save data.
function plotnsave(thlenden::ThLenDenType, params::ParamSet, nout::Int, tt::Float64)
	# Plot the shapes.
	println("\n\n\nOUTPUT NUMBER ", nout)
	nout_string = lpad(string(nout),4,string(0))
	plotfile = string(plotfolder(params),"shape",nout_string,".pdf")
	plot_curves(plotfile, thlenden.thlenvec)
	# Compute the density functions.
	getstress!(thlenden, params)
	compute_density!(thlenden, params, rotation=true)
	# Write the data to a file.
	varlabel = string("thlenden",nout_string)
	add_data(params.outfile, varlabel, thlenden)
end

# Make simple pdf plots for monitoring (not prodcution level)
function plot_curves(figname::AbstractString, thlenvec::Vector{ThetaLenType})	
	# Make figure of given height and preserve the aspect ratio.
	axlims = [1.0, 1.0]
	height = 400
	width = axlims[1]/axlims[2]*height
	# Plot the curves.
	plt = plot(xlim=(-axlims[1],axlims[1]), ylim=(-axlims[2],axlims[2]), size=(width,height),leg=false)
	for ii = 1:lastindex(thlenvec)
		thlen = thlenvec[ii]
		if thlen.len<=0
			throw("Cannot plot a curve with non-positive length.")
		end
		getxy!(thlen)
		plot!(plt, thlen.xx, thlen.yy, color="black")
	end
	savefig(plt, figname)
	return
end