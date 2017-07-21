# main.jl: The main routines to call
using Winston
include("basic.jl")
include("spectral.jl")
include("thetalen.jl")
include("RKstarter.jl")
include("plotsave.jl")
include("postprocess.jl")

#--------------- MAIN ROUTINE ---------------#
# erosion: The main routine to erode a group of bodies.
function erosion()
	# Get the input geometry, parameters, and other stuff.
	thlenden0,params,paramvec,nsteps,nout,datafolder,plotfolder = startup()
	# Begin the erosion computation with the RK starter.
	t0 = time()
	thlenden1 = RKstarter!(thlenden0, params)
	plotnsave(thlenden0,params,paramvec,datafolder,plotfolder,0.,0)
	# Enter the time loop to use the multi-step method.
	nfile = 1
	for nn = 1:nsteps
		# Compute the density function and stress for thlenden1.
		getstress!(thlenden1,params)
		# Plot and save the data in thlenden1 if appropriate.
		if mod(nn,nout)==0
			tt = nn*params.dt
			paramvec[end] = (time()-t0)/60.
			plotnsave(thlenden1,params,paramvec,datafolder,plotfolder,tt,nfile)
			nfile += 1
		end
		# Gracefully exit if all of the bodies have disappeared.
		if endof(thlenden1.thlenvec)==0; break; end
		# Advance the thlen vectors.
		advance_thetalen!(thlenden1,thlenden0,params)
	end
	# Post-process to compute the drag and other quantities.
	#postprocess("run")
	return
end

# startup: Read params and geoinfile; setup stuff.
function startup()
	# Read the parameters file.
	paramvecin = readvec("params.dat")
	# Read the input geometry file.
	geoinfile = string("../datafiles/",paramvecin[1])
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
