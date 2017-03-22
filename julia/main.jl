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
	# Set up the target points to measure u,v,p.
	ntar0 = 10; xmax = 2.8; ymax = 0.8
	ntar,xtar,ytar,utar,vtar,ptar = targets(ntar0,xmax,ymax)

	# Begin the erosion computation, time it, and save data.
	t0 = time()
	plotnsave(thlenvec0,0.,datafolder,plotfolder,0)
	# Use the Runge-Kutta starter.
	thlenvec1 = RKstarter!(thlenvec0, params)
	if nout==1
		tt = params.dt
		plotnsave(thlenvec1,tt,datafolder,plotfolder,1)
	end
	# Enter the time loop to use the multi-step method and save the data.
	nfile = 1
	for nn = 2:nsteps
		utar,vtar,ptar = stokes!(thlenvec1,params,ntar,xtar,ytar)
		advance_thetalen!(thlenvec1,thlenvec0,params)
		tt = nn*params.dt
		# Save the data.
		if mod(nn,nout)==0
			plotnsave(thlenvec1,tt,datafolder,plotfolder,nfile)
			nfile += 1
		end
		# Time the computation and write it to the params file.
		paramvec[end] = (time()-t0)/60.
		writeparams(paramsoutfile,paramvec)
	end
	return
end

# targets: Set up the target points to measure velocity and pressure: u,v,p.
function targets(nn::Integer, xmax::Float64, ymax::Float64)
	# Make the grid.
	ytar = collect(linspace(-ymax,ymax,nn))
	ytar = [ytar; ytar]
	xtar = ones(Float64,nn)
	xtar = xmax*[-xtar; xtar]
	# Initialize u,v,p at target points.
	utar,vtar,ptar = [zeros(Float64,2*nn) for ii=1:3]
	return 2*nn,xtar,ytar,utar,vtar,ptar
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

#erosion()
