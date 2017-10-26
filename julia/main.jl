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
	thlenden,params,paramvec,tfin,nout,datafolder,plotfolder = startup()
	t0 = time()
	plotnsave(thlenden,params,paramvec,datafolder,plotfolder,0.,0)
	
	# Enter the time loop to apply Runge-Kutta.
	tt = 0.; nn=0; nfile = 0;
	while(tt <= tfin && endof(thlenden.thlenvec) > 0)
		# Advance the variables forward one timestep with RK4.
		println("\n\n\n\nTIME STEP ", nn+1)
		thlenden, dt = rungekutta4(thlenden, params)
		tt += dt; nn += 1
		# Plot and save the data if appropriate.
		if mod(nn,nout)==0
			cputime = (time()-t0)/60.
			paramvec[end] = cputime
			# HERE
			# Send cputime, tt, nfile


			nfile += 1
			plotnsave(thlenden,params,paramvec,datafolder,plotfolder,tt,nfile)
		end
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
	epsilon = epsfac/npts
	sigma = sigfac/npts
	# Save params and paramvec.
	params = ParamType(dt,epsilon,sigma,nouter,ifmm,fixarea)
	paramvec = [paramvecin; cntout; 0.]
	# Create the folders for saving the data and plotting figures.
	datafolder = "../datafiles/run/" 
	newfolder(datafolder)
	plotfolder = "../figs/" 
	newfolder(plotfolder)
	return thlenden0,params,paramvec,tfin,cntout,datafolder,plotfolder
end
