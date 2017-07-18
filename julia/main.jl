# main.jl: The main routines to call
using Winston
include("basic.jl")
include("spectral.jl")
include("thetalen.jl")
include("RKstarter.jl")
include("plotsave.jl")

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
		# Advance the thlen vectors.
		advance_thetalen!(thlenden1,thlenden0,params)
		# Gracefully exit if all of the bodies have disappeared.
		if endof(thlenden1.thlenvec)==0; break; end
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

#--------------- POST PROCESSING ---------------#
# postprocess: Use the saved data to compute stuff.
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
		thlenden = new_thlenden(thlenvec,density)

		#--------------------------------------#
		# Compute velocity and pressure at a set of target points.
		xlocs = [-2.5, -2.0, -1.5, 1.5, 2.0, 2.5]
		ylocs = collect(-0.8: 0.2: 0.8)
		targets = setuptargets(xlocs,ylocs)
		compute_velpress_targets!(thlenden,targets,nouter)
		# Save the output to a data file.
		targfile = string(datafolder,"targs",cntstr,".dat")
		iostream = open(targfile, "w")
		label = string("# Stuff")
		writedlm(iostream, [label; targets.xtar; targets.ytar; 
			targets.utar; targets.vtar; targets.ptar])
		close(iostream)

		#--------------------------------------#
		# Compute the permeability
		println()
		rt1,rm1 = resistivity(thlenden, nouter, 1.5)
		rt2,rm2 = resistivity(thlenden, nouter, 2.0)
		rt3,rm3 = resistivity(thlenden, nouter, 2.5)

		#--------------------------------------#
		# Compute the drag on each body.
		# Note: the stress is not smoothed and absolute value is not taken.
		pressvec = compute_pressure(thlenden,nouter)
		tauvec = compute_stress(thlenden,nouter)
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

# setuptargets: Set up the target points.
function setuptargets(xlocs::Vector{Float64}, ylocs::Vector{Float64})
	nx = endof(xlocs)
	ny = endof(ylocs)
	ntargs = nx*ny
	xtar = ones(Float64,ntargs)
	ytar = ones(Float64,ntargs)
	for nn=1:nx
		n1,n2 = n1n2(ny,nn)
		xtar[n1:n2] = xlocs[nn]
		ytar[n1:n2] = ylocs
	end
	targets = TargetsType(evec(), evec(), evec(), evec(), evec())
	targets.xtar = xtar
	targets.ytar = ytar
	return targets
end

# resistivity: Compute the resistivity/permeability of the porous matrix.
function resistivity(thlenden::ThLenDenType, nouter::Int, x0::Float64)
	# Set up targets points on a y-grid for midpoint rule.
	nypts = 13
	dy = 2./nypts
	ylocs = collect(-1+0.5*dy: dy: 1-0.5*dy)
	# Target points for plus/minus x0.
	tarp = setuptargets([x0],ylocs)
	tarm = setuptargets([-x0],ylocs)
	# Comopute the velocities and pressures on each set of target points.
	compute_velpress_targets!(thlenden,tarp,nouter)
	compute_velpress_targets!(thlenden,tarm,nouter)
	# Compute the cross-sectional average pressure and discharge
	pplus = mean(tarp.ptar)
	pminus = mean(tarm.ptar)
	qplus = mean(tarp.utar)
	qminus = mean(tarm.utar)
	qavg = 0.5*(qplus+qminus)
	#= The discharge should be exactly the same at any location x.
	So check that it is the same at x0 and -x0. =#
	qreldiff = (qplus-qminus)/qavg
	assert(qreldiff < 1e-6)
	# Calculate the total resistivity
	rtot = (pminus - pplus)/(2*x0*qavg)
	#= Estimate the resistivity of the porous matrix region by
	an area-weighted average. =#
	Amat = 4.
	Abuf = 4*(x0-1)
	Atot = Amat+Abuf
	rbuf = 3.
	rmat = (rtot*Atot - rbuf*Abuf)/Amat
	# In the porous matrix region, subtract the contribution from the walls.
	rmat = rmat - 3.
	return rtot,rmat
end
