# main.jl: The main routines to call
using Winston
include("basic.jl")
include("spectral.jl")
include("thetalen.jl")
include("ioroutines.jl")
include("postprocess.jl")

#--------------- MAIN ROUTINE ---------------#
# erosion: The main routine to erode a group of bodies.
function erosion(dt::Float64 = -1.)
	# Get the input geometry, parameters, and other stuff.
	thlenden,params = startup()
	# Modify dt based on input to erosion if necessary.
	dt > 0? params.dt = dt : 0
	# Save the output at t=0.
	nn=0; nfile = 0; tt = 0.;
	plotnsave(nfile,tt,thlenden,params)
	# Enter the time loop to apply Runge-Kutta.
	while(tt < params.tfin - 0.1*params.dt && endof(thlenden.thlenvec) > 0)
		# Advance the variables forward one timestep with RK4.
		nn += 1
		println("\n\n\nTIME STEP ", nn)
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

# startup: Read params and geoinfile; setup stuff.
function startup()
	# Read the parameters file.
	paramvec = readvec("params.dat")
	# Read the input geometry file.
	geoinfile = string("../geometries/",paramvec[1])
	thlenvec0 = read_thlen_file(geoinfile)
	thlenden0 = new_thlenden(thlenvec0)
	npts,nbods = getnvals(thlenvec0)
	# Define the object params.
	params = getparams(paramvec,npts)
	# Create new data folders
	datafolder, plotfolder = getfoldernames()
	newfolder(datafolder)
	newfolder(plotfolder)
	return thlenden0,params
end
# function getparams: Define the object of parameters.
function getparams(paramvec::Vector, npts::Int)
	# Read the parameters and calculate needed quantities.
	epsfac,sigfac,dt,dtout,tfin,nouter,ifmm,fixarea = paramvec[2:9]
	epsilon = epsfac/npts
	sigma = sigfac/npts
	cntout = max(round(Int,dtout/dt),1)
	cput0 = time()
	# Save the parameters in an object.
	params = ParamType(dt,epsilon,sigma,nouter,ifmm,fixarea,npts,tfin,cntout,cput0)
	return params
end
