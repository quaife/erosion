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
	thlenden,params = startup()
	# Save the output at t=0.
	nn=0; nfile = 0; tt = 0.;
	plotnsave(nfile,tt,thlenden,params)
	
	# Enter the time loop to apply Runge-Kutta.
	while(tt <= params.tfin && endof(thlenden.thlenvec) > 0)
		# Advance the variables forward one timestep with RK4.
		nn += 1
		println("\n\n\n\nTIME STEP ", nn)
		thlenden, dt = rungekutta4(thlenden, params)
		tt += dt
		# Plot and save the data if appropriate.
		if mod(nn, params.cntout)==0
			nfile += 1
			plotnsave(nfile,tt,thlenden,params)
		end
	end
	return
end

# startup: Read params and geoinfile; setup stuff.
function startup()
	# Read the parameters file.
	paramvec = readvec("params.dat")
	# Read the input geometry file.
	geoinfile = string("../geometries/",paramvec[1])
	thlenvec0 = read_thlen_file(geoinfile)
	thlenden0 = new_thlenden(thlenvec0)
	# Read the other parameters and calculate needed quantities.
	nouter,tfin,dtout,dtfac,epsfac,sigfac,ifmm,fixarea = paramvec[2:9]
	npts,nbods = getnvals(thlenden0.thlenvec)
	dt = dtfac/npts
	cntout = round(Int,dtout/dt)
	cntout = max(cntout,1)
	epsilon = epsfac/npts
	sigma = sigfac/npts
	cput0 = time()
	# Save the parameters in an object.
	params = ParamType(dt,epsilon,sigma,nouter,ifmm,fixarea,tfin,cntout,cput0)
	# Create new data folders
	datafolder, plotfolder = getfoldernames()
	newfolder(datafolder)
	newfolder(plotfolder)
	return thlenden0,params
end
