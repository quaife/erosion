# main.jl: The main routines to call
using Winston
include("basic.jl")
include("spectral.jl")
include("thetalen.jl")
include("RKstarter.jl")
include("plotsave.jl")

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
		getstress!(thlenden1,params)
		advance_thetalen!(thlenden1,thlenden0,params)
		# Gracefully exit if all of the bodies have disappeared.
		if endof(thlenden1.thlenvec)==0; break; end
		# Plot and save the data when appropriate.
		# Note: Save thlenden0 because thlenden1 does not yet have the density-function computed.
		if mod(nn,nout)==0
			tt = nn*params.dt
			paramvec[end] = (time()-t0)/60.
			plotnsave(thlenden0,params,paramvec,datafolder,plotfolder,tt,nfile)
			nfile += 1
		end
	end
	# Post-process to compute the drag and other quantities.
	postprocess("run")
	return
end

# startup: Read params and geoinfile; setup stuff.
function startup()
	# Read the parameters file.
	paramvecin = readvec("params.dat")
	# Read the input geometry file.
	geoinfile = string("../datafiles/",paramvecin[1])
	thlenvec0 = read_thlen_file(geoinfile)
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

# postprocess:
function postprocess(foldername::AbstractString)
	# Define the data folder.
	datafolder = string("../datafiles/",foldername,"/")
	# Read the params data file.
	paramsfile = string(datafolder,"params.dat")
	paramvec = readvec(paramsfile)
	nouter = paramvec[2]
	ntimes = paramvec[end]
	# Read the data at each time step.
	for cnt=0:ntimes
		# Get the file name at each time.
		cntstr = lpad(cnt,4,0)
		geomfile = string(datafolder,"geom",cntstr,".dat")
		densityfile = string(datafolder,"density",cntstr,".dat")
		# Extract thlenvec and density.
		tt,thlenvec = read_geom_file(geomfile)
		density = readvec(densityfile)
		thlenden = ThLenDenType(thlenvec,density)
		# Compute the pressure and stress on all bodies.
		# Note: the stress is not smoothed and absolute value is not taken.
		pressvec = computepressure(thlenden,nouter)
		tauvec = computestress(thlenden,nouter)
		#--------------------------------------#
		# Compute the drag on each body.
		npts,nbods = getnvals(thlenvec)
		dragxvec,dragyvec = [zeros(Float64,nbods) for ii=1:2]
		for nn=1:nbods
			# Get the pressure and stress and boddy nn.
			n1,n2 = n1n2(npts,nn)
			press = pressvec[n1:n2]
			tau = tauvec[n1:n2]
			# Get the tangent/normal vectors and arc length increment.
			sx,sy,nx,ny = getns(thlenvec[nn].theta)
			ds = thlenvec[nn].len / npts
			# Compute the drag force.
			# Note: I believe both should be plus signs due to the conventions.
			dragxvec[nn] = sum(press.*nx + tau.*sx)*ds
			dragyvec[nn] = sum(press.*ny + tau.*sy)*ds
		end
		# Save the output to a data file.
		dragfile = string(datafolder,"drag",cntstr,".dat")
		iostream = open(dragfile, "w")
		label = string("# Drag data for ",nbods," bodies. All xdrags then all ydrags.")
		writedlm(iostream, [label; dragxvec; dragyvec])
		close(iostream)
	end
	return
end

