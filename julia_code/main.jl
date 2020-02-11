# main.jl: The main routines to call
using Plots
using DelimitedFiles
using Statistics
using FFTW
include("basic.jl")
include("makegeos.jl")
include("spectral.jl")
include("thetalen.jl")
include("ioroutines.jl")
include("postprocess.jl")


#--------------- MAIN ROUTINE ---------------#
# Dispatch to call the main routine with dt and tfin set by params file.
function erosion(paramsfile::AbstractString = "params")
	thlenden,params = startup(paramsfile)
	erosion(thlenden,params)
	postprocess(string("run_",paramsfile))
end
# Dispatch to call the main routine with dt and tfin set by the caller.
function erosion(paramsfile::AbstractString, dt::Float64, tfin::Float64)
	thlenden,params = startup(paramsfile)
	params.dt = dt; params.tfin = tfin
	erosion(thlenden,params)
end
# erosion: The main routine to erode a group of bodies.
function erosion(thlenden::ThLenDenType, params::ParamType)
	println("Running erosion with dt = ", round(params.dt, sigdigits=3), 
		" and tfin = ", round(params.tfin, sigdigits=3))
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
	cputime = round( (time()-params.cput0)/3600. , sigdigits=2)
	println("\n\n\nCOMPLETED SIMULATION")
	println("cpu time = ", cputime, " hours.\n\n")
	return thlenden,params,tt,cputime
end

# startup: Read params and geoinfile; setup stuff.
function startup(paramsfile::AbstractString)
	# Read the input geometry file.
	paramvec = readvec(string(paramsfile,".txt"))
	circfile = string("../", paramvec[1])
	npts = paramvec[2]
	thlenfile = string("../thlen_tmp.txt")
	save_thlen(circfile, thlenfile, npts)
	thlenvec0 = read_thlen_file(thlenfile)
	thlenden0 = new_thlenden(thlenvec0)
	npts1,nbods = getnvals(thlenvec0)
	assert(npts == npts1)
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
# function getparams: Define the object of parameters.
function getparams(paramsfile::AbstractString)
	# Read the parameters.
	paramvec = readvec(string(paramsfile,".txt"))
	circfile = paramvec[1]
	npts = paramvec[2]
	ibary, ifmm = Int(paramvec[3]), Int(paramvec[4])
	epsfac, sigfac, dt, dtout = paramvec[5:8]
	fixpdrop, fixarea = Bool(paramvec[9]), Bool(paramvec[10])
	tfin = paramvec[11]
	maxl, nouter = Int(paramvec[12]), Int(paramvec[13])
	# Calculate the needed quantities.
	epsilon = epsfac/npts
	sigma = sigfac/npts
	cntout = max(round(Int,dtout/dt),1)
	cput0 = time()
	# Save the parameters in an object.
	params = ParamType(dt,epsilon,sigma,nouter,ifmm,ibary,maxl,
		fixarea,fixpdrop,npts,tfin,cntout,cput0,geofile,paramsfile)
	return params
end
