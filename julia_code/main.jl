# main.jl: The main routines to simulate erosion
# Saves the output data in a jld2 file.

using JLD2
using Plots

include("basic.jl")
include("callFortran.jl")
include("spectral.jl")
include("thetalen.jl")
 
 #--------------- TINY ROUTINES ---------------#
 # Set the plot folder.
plotfolder(params::ParamSet) = string("../zFigs-",params.label,"/")

# Add a variable incrementally to a jld data file.
function add_data(filename::AbstractString, varlabel::AbstractString, var)
	jldopen(filename, "r+") do file
		write(file, varlabel, var)
	end
end

# Convert an integer to a string with zero-padding.
nstr(nn::Int) = lpad(string(nn), 4, string(0))
# The data label for thlenden
thlabel(nn) = string("thlenden", nstr(nn))
#-------------------------------------------------#

#--------------- SMALL ROUTINES ---------------#
# Convert the circle data to thlen data.
function circ2thlen(npts::Int, rad::Float64, xc::Float64, yc::Float64)
	alpha = getalpha(npts)
	theta = 0.5*pi .+ alpha
	len = 2*pi*rad
	return ThetaLenType(theta,len,xc,yc,NaN)
end
# circs2thlenden: Initialize thlenden from the input circle file.
function circs2thlenden(params::ParamSet)
	circdata = readvec(params.infile)
	nbods = round(Int, popfirst!(circdata))
	thlenvec = Array{ThetaLenType}(undef, 0)
	for nn = 1:nbods
		rad, xc, yc = [popfirst!(circdata) for i=1:3]
		thlen = circ2thlen(params.npts, rad, xc, yc)
		push!(thlenvec, thlen)
	end
	return new_thlenden(thlenvec)
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
		xx,yy = getxy(thlen)
		plot!(plt, xx, yy, color="black")
	end
	savefig(plt, figname)
	return
end

# Plot and save data.
function plotnsave(thlenden::ThLenDenType, params::ParamSet, nout::Int)
	# Plot the shapes.
	println("\n\n\nOUTPUT NUMBER ", nout)
	plotfile = string(plotfolder(params),"shape",nstr(nout),".pdf")
	plot_curves(plotfile, thlenden.thlenvec)
	# Compute the density functions.
	compute_density!(thlenden, params)
	compute_density!(thlenden, params, rotation=true)
	# Add the thlenden data to the data file.
	add_data(params.outfile, thlabel(nout), thlenden)
end
#-------------------------------------------------#


#--------------- MAIN ROUTINES ---------------#
# The routine to erode a group of bodies.
function erosion(params::ParamSet)
	# Initialize.
	thlenden = circs2thlenden(params)
	# If the plot folder exists, delete it, then create a new folder.
	pfolder = plotfolder(params)
	if isdir(pfolder) rm(pfolder; recursive=true) end
	mkdir(pfolder)
	nn, nout = 0, 0
	# Enter the time loop to apply Runge-Kutta.
	while(thlenden.tt < params.tfin && length(thlenden.thlenvec) > 0)
		# Plot and save the data if appropriate.
		if mod(nn, params.outstride) == 0
			plotnsave(thlenden, params, nout)
			nout += 1
		end
		# Advance the variables forward one timestep with RK4.
		nn += 1
		println("\n\nTIME STEP ", nn)
		thlenden = rungekutta2(thlenden, params)
	end
	# Plot and save one last time with zero bodies.
	plotnsave(thlenden, params, nout)
	add_data(params.outfile, "noutputs", nout)
end

# The main routine to call erosion.
function main(params::ParamSet)
	# Initialize the output jld2 file and save the parameters.
	jldsave(params.outfile; params)
	# Run the erosion simulation.
	println("\nBEGINNING EROSION SIMULATION")
	cputime = @elapsed	erosion(params)
	# Save the CPU time of the simulation.
	cpu_hours = round(cputime/3600, sigdigits=3)
	add_data(params.outfile, "cpu_hours", cpu_hours)
	println("\n\n\nCOMPLETED EROSION SIMULATION")
	println("cpu time = ", cpu_hours, " hours.\n\n")
end
#-------------------------------------------------#