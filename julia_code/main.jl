# main.jl: The main routines to call
using Plots
using DelimitedFiles
using Statistics
using FFTW
include("basic.jl")
##include("makegeos.jl")
include("spectral.jl")
include("thetalen.jl")
include("ioroutines.jl")


sig(var::AbstractFloat, sigdig::Int=3) = round(var,sigdigits=sigdig)


#--------------- MAIN ROUTINE ---------------#
# The main routine to erode a group of bodies.
function erosion(params::paramset)


	thlenden = get_thlenden(params)



	println("Running erosion: dt = ",round(params.dt),"; tfin = ", sig(params.tfin))
	# Save the output at t=0.
	nn=0; nfile = 0; tt = 0.;
	plotnsave(nfile,tt,thlenden,params)
	# Enter the time loop to apply Runge-Kutta.
	while(tt < params.tfin - 0.1*params.dt && length(thlenden.thlenvec) > 0)
		# Print statements
		nn += 1
		println("\n\n\nTIME STEP ", nn)
		println("t/tfin = ", round(tt/params.tfin, sigdigits=3))
		cpumins = round( (time()-params.cput0)/60. , sigdigits=2)
		println("cpu time = ", cpumins, " minutes.")
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

	return thlenden,params,tt,cputime
end


# Get stuff.
function get_thlenden(params::paramset)
	# Read the input circle file.
	npts = params.npts
	circdata = readvec(params.infile)
	nbods = round(Int, circdata[1])
	deleteat!(circdata,1)

	for nn = 1:nbods
		rad, xc, yc = circdata[1], circdata[2], circdata[3]
		thlen = circthlen(npts, rad, xc, yc)
		deleteat!(circdata,1:3)
	end


	save_thlen(circfile, thlenfile, npts)
	thlenvec0 = read_thlen_file(thlenfile)
	thlenden0 = new_thlenden(thlenvec0)


end









# startup: Read params and geoinfile; setup stuff.
function startup(paramsfile::AbstractString)
	# Read the input geometry file.
	paramvec = readvec(string(paramsfile,".txt"))
	circfile = string("../", paramvec[1])
	npts = paramvec[2]
	thlenfile = string("../input_geos/thlen_tmp.txt")
	save_thlen(circfile, thlenfile, npts)
	thlenvec0 = read_thlen_file(thlenfile)
	thlenden0 = new_thlenden(thlenvec0)
	npts1,nbods = getnvals(thlenvec0)
	@assert(npts == npts1)
	# Define the object params.
	params = getparams(paramsfile)
	# Create new data folders
	datafolder, plotfolder = getfoldernames(params.paramsfile)
	newfolder(datafolder)
	newfolder(plotfolder)
	# Save the input parameters file.
	saveparamsfile = string(datafolder,"aparams.txt")
	writedata(paramvec,saveparamsfile)
	return thlenden0,params
end
