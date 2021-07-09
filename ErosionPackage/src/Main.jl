# MAIN GOAL: The main routines to simulate erosion.
# The output is saved in a jld2 file.
# Convention: nn indexes the timestep; bod = 1:nbods indexes the bodies.

using JLD2
using Plots
using DelimitedFiles

using .ThetaLen: ParamSet

#-------------------------------------------------#

# USED IN ONLY ONE PLACE, THINK ABOUT IT...
# Read a vector from a text file.
function readvec(filename::AbstractString)
	iostream = open(filename, "r")
	invec = readdlm(iostream, comments=true)[:,1]
	close(iostream)
	return invec
end

#--------------- TINY ROUTINES ---------------#
# Set files and folders: the input data file, the output temporary data file,
# the final output file, and the output plot folder.
infile(params::ParamSet) = string(params.infolder, params.label, ".circ")
tempfile(params::ParamSet) = string("temp_data-", params.label, ".jld2")
outfile(params::ParamSet) = string("raw_data-", params.label, ".jld2")
plotfolder(params::ParamSet) = string("zFigs-", params.label, "/")

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
# Initialize thlenden from the input circle file.
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
function plot_curves(thlenvec::Vector{ThetaLenType}, params::ParamSet, nout::Int)
	println("\n\n\nOUTPUT NUMBER ", nout)
	plotfile = string(plotfolder(params), "shape", nstr(nout), ".pdf")	
	# Make figure of given height and preserve the aspect ratio.
	axlims = [1.0, 1.0]
	height = 400
	width = axlims[1]/axlims[2]*height
	plt = plot(xlim=(-axlims[1],axlims[1]), ylim=(-axlims[2],axlims[2]), 
		size=(width,height), leg=false)
	# Plot the curves.
	for bod = 1:length(thlenvec)
		thlen = thlenvec[bod]
		if thlen.len<=0; throw("Cannot plot a curve with non-positive length."); end
		xx, yy = getxy(thlen)
		plot!(plt, xx, yy, color="black")
	end
	savefig(plt, plotfile)
	return
end

# Compute both density functions, regular and rotated, and add thlenden to thldvec
function update_thldvec!(thldvec::Vector{ThLenDenType}, thlenden::ThLenDenType, params::ParamSet)
	# Compute the density functions, both regular and rotated.
	compute_density!(thlenden, params)
	compute_density!(thlenden, params, rotation=true)
	push!(thldvec, thlenden)
	return
end
#-------------------------------------------------#


#--------------- MAIN ROUTINES ---------------#
# The routine to erode a group of bodies.
# This routine returns thldvec, the vector containing thlenden at all time values.
function erosion_sim(params::ParamSet)
	# Initialize the temporary output file by saving the parameters.
	jldsave(tempfile(params); params)
	
	# Initialize the geometry.
	thlenden = circs2thlenden(params)
	
	# Initialize other stuff.
	println("\nBEGINNING EROSION SIMULATION")
	pfolder = plotfolder(params)
	if isdir(pfolder) rm(pfolder; recursive=true) end; mkdir(pfolder)
	thldvec = Vector{ThLenDenType}(undef, 0)
	nn, nout = 0, 0

	# Enter the time loop to apply Runge-Kutta.
	while(thlenden.tt < params.tfin && length(thlenden.thlenvec) > 0)
		
		# Save and plot the data if appropriate.
		if mod(nn, params.outstride) == 0
			# Important: Update the vector of saved thlenden values		.	
			update_thldvec!(thldvec, thlenden, params)
			# Also plot the curves for monitoring.
			plot_curves(thlenden.thlenvec, params, nout)
			# Also incrementally add thlenden to the temporary data file in case of crash.
			add_data(tempfile(params), thlabel(nout), thlenden)
			nout += 1
		end
		# Advance the variables forward one timestep with RK4.
		nn += 1
		println("\n\nTIME STEP ", nn)
		thlenden = rungekutta2(thlenden, params)
	end

	# Update one last time with zero bodies and return thldvec.
	update_thldvec!(thldvec, thlenden, params)
	println("\n\n\nCOMPLETED EROSION SIMULATION")
	return thldvec
end

# The main routine to call the erosion simulation.
#= Note: This routine was created separately primarily so that 
the simulation erosion() routine could be timed. =#
function run_erosion(params::ParamSet)
	cputime = @elapsed	thldvec = erosion_sim(params)

	# Calculate the CPU time of the simulation and print it.
	cpu_hours = round(cputime/3600, sigdigits=3)
	println("cpu time = ", cpu_hours, " hours.\n\n")

	# Save the data to the main output file.
	jldsave(outfile(params); params, thldvec, cpu_hours)
end