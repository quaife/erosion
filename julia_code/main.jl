# MAIN GOAL: The main routines to simulate erosion.
# The output is saved in a jld2 file.
# Convention: nn indexes the timestep, bod = 1:nbods indexes the bodies.


using TimeStep


using JLD2
using Plots

# USED IN ONLY ONE PLACE, THINK ABOUT IT...
# basic.jl: Basic routines such as datatypes.
using DelimitedFiles
# readvec: Read a vector from a text file.
function readvec(filename::AbstractString)
	iostream = open(filename, "r")
	invec = readdlm(iostream, comments=true)[:,1]
	close(iostream)
	return invec
end





 #--------------- TINY ROUTINES ---------------#
 # Set files and folders: the input data file, the output temporary data file,
 # the final output file, and the output plot folder.
infile(params::ParamSet) = string("../", params.infolder, params.label, ".circ")
tempfile(params::ParamSet) = string("../temp_data-", params.label, ".jld2")
outfile(params::ParamSet) = string("../raw_data-", params.label, ".jld2")
plotfolder(params::ParamSet) = string("../zFigs-", params.label, "/")

# Add a variable incrementally to a jld data file.
function add_data(filename::AbstractString, varlabel::AbstractString, var)
	jldopen(filename, "r+") do file
		write(file, varlabel, var)
	end
end

# Convert an integer to a string with zero-padding.
nstr(nn::Int) = lpad(string(nn), 4, string(0))
# The data label for thlenden.
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
	circdata = readvec(infile(params))
	nbods = round(Int, popfirst!(circdata))
	thlenvec = Array{ThetaLenType}(undef, 0)
	for bod = 1:nbods
		rad, xc, yc = [popfirst!(circdata) for i=1:3]
		thlen = circ2thlen(params.npts, rad, xc, yc)
		push!(thlenvec, thlen)
	end
	return new_thlenden(thlenvec)
end

# Make simple pdf plots for monitoring (not production level).
function plot_curves(figname::AbstractString, thlenvec::Vector{ThetaLenType})	
	# Make figure of given height and preserve the aspect ratio.
	axlims = [1.0, 1.0]
	height = 400
	width = axlims[1]/axlims[2]*height
	plt = plot(xlim=(-axlims[1],axlims[1]), ylim=(-axlims[2],axlims[2]), 
		size=(width,height), leg=false)
	# Plot the curves.
	for ii = 1:lastindex(thlenvec)
		thlen = thlenvec[ii]
		if thlen.len<=0
			throw("Cannot plot a curve with non-positive length.")
		end
		xx, yy = getxy(thlen)
		plot!(plt, xx, yy, color="black")
	end
	savefig(plt, figname)
	return
end

# Plot the data and incrementally save it to the temporary output file.
function plotnsave(thlenden::ThLenDenType, params::ParamSet, nout::Int)
	# Plot the shapes.
	println("\n\n\nOUTPUT NUMBER ", nout)
	plotfile = string(plotfolder(params),"shape",nstr(nout),".pdf")
	plot_curves(plotfile, thlenden.thlenvec)
	# Compute the density functions.
	compute_density!(thlenden, params)
	compute_density!(thlenden, params, rotation=true)
	# Add the thlenden data to the data file.
	add_data(tempfile(params), thlabel(nout), thlenden)
end
#-------------------------------------------------#


#--------------- MAIN ROUTINES ---------------#
# The routine to erode a group of bodies.
# This routine saves the output in the variable thldvec.
function erosion!(params::ParamSet, thldvec::Vector{ThLenDenType})
	# Initialize the temporary output file by saving the parameters.
	jldsave(tempfile(params); params)

	# Initialize the geometry and the plot folder.
	thlenden = circs2thlenden(params)
	pfolder = plotfolder(params)
	if isdir(pfolder) rm(pfolder; recursive=true) end; mkdir(pfolder)
	
	# Enter the time loop to apply Runge-Kutta.
	nn, nout = 0, 0
	while(thlenden.tt < params.tfin && length(thlenden.thlenvec) > 0)
		# Plot and save the data if appropriate.
		if mod(nn, params.outstride) == 0			
			plotnsave(thlenden, params, nout)
			push!(thldvec, thlenden)
			nout += 1
		end
		# Advance the variables forward one timestep with RK4.
		nn += 1
		println("\n\nTIME STEP ", nn)
		thlenden = rungekutta2(thlenden, params)
	end
	# Plot and save one last time with zero bodies.
	plotnsave(thlenden, params, nout)
	push!(thldvec, thlenden)
end


# The main routine to call erosion.
#= Note: This routine was created separately primarily so that 
the erosion() routine could be timed. =#
function main(params::ParamSet)
	# Run the erosion simulation and time it.
	println("\nBEGINNING EROSION SIMULATION")
	thldvec = Vector{ThLenDenType}(undef, 0)
	cputime = @elapsed	erosion!(params, thldvec)

	# Calculate the CPU time of the simulation and print it.
	cpu_hours = round(cputime/3600, sigdigits=3)
	println("\n\n\nCOMPLETED EROSION SIMULATION")
	println("cpu time = ", cpu_hours, " hours.\n\n")

	# Save the data to the main output file.
	jldsave(outfile(params); params, thldvec, cpu_hours)
end
#-------------------------------------------------#