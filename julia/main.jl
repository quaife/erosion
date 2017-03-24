# main.jl: The main routines to call
#################### Includes ####################
using Winston
include("basic.jl")
include("spectral.jl")
include("thetalen.jl")
include("misc.jl")
##################################################

# erosion: The main routine to erode a group of bodies.
function erosion()
	# Get the input geometry and parameters.
	thlenvec0,params,paramvec,nsteps,nout = getparams()
	# Create the folders for saving the data and plotting figures
	datafolder = "../datafiles/run/"; newfolder(datafolder)
	plotfolder = "../figs/"; newfolder(plotfolder)
	paramsoutfile = string(datafolder,"params.dat")

	# Begin the erosion computation with the RK starter.
	t0 = time()
	plotnsave(thlenvec0,0.,datafolder,plotfolder,0)
	thlenvec1 = RKstarter!(thlenvec0, params)
	if nout==1; plotnsave(thlenvec1,params.dt,datafolder,plotfolder,1); end
	
	# Enter the time loop to use the multi-step method.
	nfile = 1
	for nn = 2:nsteps
		

		stokes!(thlenvec1,params)
		advance_thetalen!(thlenvec1,thlenvec0,params)
		
		# Plot and save the data when appropriate.
		if mod(nn,nout)==0
			# Plot and save the data.
			tt = nn*params.dt
			plotnsave(thlenvec1,tt,datafolder,plotfolder,nfile)
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
	npts = endof(thlenvec0[1].theta)
	# Read the rest of the parameters.
	tfin,dtout,dtfac,epsfac,sigfac,ifmm,fixarea = paramvecin[2:8]
	# Calculate the needed parameters.
	dt = dtfac/npts
	cntout = round(Int,dtout/dt)
	cntout = max(cntout,1)
	nsteps = round(Int,tfin/dt)
	epsilon = epsfac/npts
	sigma = sigfac/npts
	dtoutexact = cntout*dt
	# Save params and paramvec
	params = ParamType(dt,epsilon,sigma,ifmm,fixarea)
	paramvec = [paramvecin; dtoutexact; cntout; 0.]
	return thlenvec0,params,paramvec,nsteps,cntout
end
