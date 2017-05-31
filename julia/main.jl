# main.jl: The main routines to call
#################### Includes ####################
using Winston
include("basic.jl")
include("spectral.jl")
include("thetalen.jl")
include("RKstarter.jl")
include("misc.jl")
##################################################

# erosion: The main routine to erode a group of bodies.
function erosion()
	# Get the input geometry and parameters.
	thlenden0,params,paramvec,nsteps,nout = getparams()
	dt = params.dt

	# Create the folders for saving the data and plotting figures
	datafolder = "../datafiles/run/"; newfolder(datafolder)
	plotfolder = "../figs/"; newfolder(plotfolder)
	paramsoutfile = string(datafolder,"params.dat")
	# Begin the erosion computation with the RK starter.
	t0 = time()
	thlenden1 = RKstarter!(thlenden0, params)
	plotnsave(thlenden0,params,datafolder,plotfolder,0.,0)
	# Plot and save the data if appropriate.
	nfile = 1
	if nout==1
		plotnsave(thlenden1,params,datafolder,plotfolder,dt,1)
		nfile += 1
	end
	# Enter the time loop to use the multi-step method.
	for nn = 2:nsteps
		getstress!(thlenden1,params)
		advance_thetalen!(thlenden1,thlenden0,params)
		# Plot and save the data when appropriate.
		if mod(nn,nout)==0
			# Plot and save the data.
			tt = nn*dt
			plotnsave(thlenden1,params,datafolder,plotfolder,tt,nfile)
			# Time the computation and write it to the params file.
			paramvec[end] = (time()-t0)/60.
			writeparams(paramsoutfile,paramvec)
			nfile += 1
		end
	end
	return
end

function getparams()
	# Read the parameters file.
	iostream = open("params.dat", "r")
	paramvecin = readdlm(iostream)
	close(iostream)
	# Read the input geometry file.
	geoinfile = string(paramvecin[1])
	thlenvec0 = readthlenfile(string("../datafiles/",geoinfile))
	thlenden0 = ThLenDenType(thlenvec0,evec())
	# Read the other parameters and calculate needed quantities.
	nouter,tfin,dtout,dtfac,epsfac,sigfac,ifmm,fixarea = paramvecin[2:9]
	npts = endof(thlenvec0[1].theta)
	dt = dtfac/npts
	cntout = round(Int,dtout/dt)
	cntout = max(cntout,1)
	nsteps = round(Int,tfin/dt)
	epsilon = epsfac/npts
	sigma = sigfac/npts
	dtoutexact = cntout*dt
	# Save params and paramvec
	params = ParamType(dt,epsilon,sigma,nouter,ifmm,fixarea)
	paramvec = [paramvecin; dtoutexact; cntout; 0.]
	return thlenden0,params,paramvec,nsteps,cntout
end
