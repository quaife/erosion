# main.jl: The main routines to call
using Winston
include("basic.jl")
include("spectral.jl")
include("thetalen.jl")
include("ioroutines.jl")
include("postprocess.jl")

#--------------- MAIN ROUTINE ---------------#
# erosion: The main routine to erode a group of bodies.
function erosion()
	# Get the input geometry, parameters, and other stuff.
	thlenden,params,paramvec,nsteps,nout,datafolder,plotfolder = startup()
	# Begin the erosion computation with the RK starter.
	t0 = time()
	plotnsave(thlenden,params,paramvec,datafolder,plotfolder,0.,0)
	
	# Enter the time loop to apply Runge-Kutta.
	nfile = 1; tt = 0.
	for nn = 1:nsteps
		# Advance the variables forward one timestep with RK4.
		println("\n\n\n\nTIME STEP ", nn)
		thlenden, dt = rungekutta4(thlenden, params)
		tt += dt
		# Plot and save the data if appropriate.
		if mod(nn,nout)==0
			paramvec[end] = (time()-t0)/60.
			plotnsave(thlenden,params,paramvec,datafolder,plotfolder,tt,nfile)
			nfile += 1
		end
		# Gracefully exit if all of the bodies have disappeared.
		if endof(thlenden.thlenvec)==0; break; end
	end
	return
end

# startup: Read params and geoinfile; setup stuff.
function startup()
	# Read the parameters file.
	paramvecin = readvec("params.dat")
	# Read the input geometry file.
	geoinfile = string("../geometries/",paramvecin[1])
	thlenvec0 = read_thlen_file(geoinfile)
	thlenden0 = new_thlenden(thlenvec0)
	# Read the other parameters and calculate needed quantities.
	nouter,tfin,dtout,dtfac,epsfac,sigfac,ifmm,fixarea = paramvecin[2:9]
	npts,nbods = getnvals(thlenden0.thlenvec)
	dt = dtfac/npts
	cntout = round(Int,dtout/dt)
	cntout = max(cntout,1)
	nsteps = round(Int,tfin/dt)
	epsilon = epsfac/npts
	sigma = sigfac/npts
	dtoutexact = cntout*dt
	# Save params and paramvec.
	params = ParamType(dt,epsilon,sigma,nouter,ifmm,fixarea)
	paramvec = [paramvecin; dtoutexact; cntout; 0.]
	# Create the folders for saving the data and plotting figures.
	datafolder = "../datafiles/run/" 
	newfolder(datafolder)
	plotfolder = "../figs/" 
	newfolder(plotfolder)
	return thlenden0,params,paramvec,nsteps,cntout,datafolder,plotfolder
end
