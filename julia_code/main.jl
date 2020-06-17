# main.jl: The main routines to call
using JLD
using DelimitedFiles
using Plots
using Statistics
using FFTW
include("basic.jl")
include("spectral.jl")
include("thetalen.jl")
 


# Shorthand to round a float to a number of significant digits.
sig(var::AbstractFloat, sigdig::Int=3) = round(var,sigdigits=sigdig)

# Add a variable incrementally to a jld data file.
function add_data(filename::AbstractString, varlabel::AbstractString, var)
	jldopen(filename, "r+") do file
	file["varlabel"] = var
	end
end

# Set the plot folder.
plotfolder(params::ParamSet) = string("../zFigs-",params.label,"/")
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
function main(params::ParamSet)
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
		println("\n\n\nTIME STEP ", nn)
		thlenden, dt = rungekutta2(thlenden, params)
		tt += dt
	end
	# Plot and save one last time with zero bodies.
	plotnsave(thlenden,params,nout,tt)
end





# Get initial thlenden.
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
function plotnsave(thlenden::ThLenDenType, params::ParamSet, nfile::Int, tt::Float64)
	# Plot the shapes.
	println("\n\n\nOUTPUT NUMBER ", nfile)
	nfilestr = lpad(string(nfile),4,string(0))
	plotfile = string(plotfolder(params),"shape",nfilestr,".pdf")
	plot_curves(plotfile, thlenden.thlenvec)
	# Compute the density functions.
	getstress!(thlenden, params)
	compute_density!(thlenden, params, rotation=true)



	# Write the data to a file.
	add_data(params.outfile, varlabel, thlenden)


	# TO DO: Need to push! the new thlenden to the thlendenvec

	## save_geo_density(tt,thlenden,geomfile,densityfile)
	## save_pinfo(params,nfile,pinfofile)

	return
end



# plotcurve: Plot multiple curves from the theta-len values.
function plot_curves(figname::AbstractString, thlenvec::Vector{ThetaLenType})	
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




