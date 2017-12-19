# main.jl: The main routines to call
using Winston
include("basic.jl")
include("spectral.jl")
include("thetalen.jl")
include("ioroutines.jl")
include("postprocess.jl")


#--------------- MAIN ROUTINE ---------------#
# Dispatch to call the main routine with dt and tfin set by params file.
function erosion(paramsfile::AbstractString = "params.in")
	thlenden,params = startup(paramsfile)
	erosion(thlenden,params)
end
# Dispatch to call the main routine with dt and tfin set by the caller.
function erosion(paramsfile::AbstractString, dt::Float64, tfin::Float64)
	thlenden,params = startup(paramsfile)
	params.dt = dt; params.tfin = tfin
	erosion(thlenden,params)
end
# erosion: The main routine to erode a group of bodies.
function erosion(thlenden::ThLenDenType, params::ParamType)
	println("Running erosion with dt = ", signif(params.dt,3), 
		" and tfin = ", signif(params.tfin,3))
	# Save the output at t=0.
	nn=0; nfile = 0; tt = 0.;
	plotnsave(nfile,tt,thlenden,params)
	# Enter the time loop to apply Runge-Kutta.
	while(tt < params.tfin - 0.1*params.dt && endof(thlenden.thlenvec) > 0)
		# Advance the variables forward one timestep with RK4.
		nn += 1
		println("\n\n\nTIME STEP ", nn)
		println("t/tfin = ", signif(tt/params.tfin,3))
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
	cputime = round( (time()-params.cput0)/60. , 2)
	println("\n\n\nCOMPLETED SIMULATION")
	println("cpu time = ", cputime, " minutes.\n\n")
	return thlenden,params,tt
end

# startup: Read params and geoinfile; setup stuff.
function startup(paramsfile::AbstractString)
	# Read the input geometry file.
	paramvec = readvec(paramsfile)
	geoinfile = string("../geometries/",paramvec[1],".in")
	thlenvec0 = read_thlen_file(geoinfile)
	thlenden0 = new_thlenden(thlenvec0)
	npts,nbods = getnvals(thlenvec0)
	# Define the object params.
	params = getparams(paramsfile,npts)
	# Create new data folders
	datafolder, plotfolder = getfoldernames(params.paramsfile)
	newfolder(datafolder)
	newfolder(plotfolder)
	# Save the input parameters file.
	saveparamsfile = string(datafolder,"aparams.in")
	writedata(paramvec,saveparamsfile)
	return thlenden0,params
end
# function getparams: Define the object of parameters.
function getparams(paramsfile::AbstractString, npts::Int)
	# Read the parameters and calculate needed quantities.
	paramvec = readvec(paramsfile)
	geofile,epsfac,sigfac,dt,dtout,tfin,nouter,ifmm = paramvec[1:8]
	fixarea,fixpdrop = Bool(paramvec[9]),Bool(paramvec[10])
	epsilon = epsfac/npts
	sigma = sigfac/npts
	cntout = max(round(Int,dtout/dt),1)
	cput0 = time()
	# Save the parameters in an object.
	params = ParamType(dt,epsilon,sigma,nouter,ifmm,
		fixarea,fixpdrop,npts,tfin,cntout,cput0,geofile,paramsfile)
	return params
end
