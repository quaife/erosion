# main.jl: The main routines to call
#################### Includes ####################
using Winston
include("basic.jl")
include("spectral.jl")
include("thetalen.jl")
include("misc.jl")
##################################################

# erosion: The main routine to erode a group of bodies.
function erosion(thleninput::AbstractString)
	# Read the input geometry file in thetlen form.
	thlenvec0 = readthlenfile(string("../datafiles/",thleninput))
	npts = endof(thlenvec0[1].theta)
	# Read the parameters from the input data file.
	invec = readparams()
	tfin = invec[1]
	dtout = invec[2]
	dtfac = invec[3]
	epsfac = invec[4]
	sigfac = invec[5]
	lenevo = invec[6]
	# Calculate the needed parameters.
	dt = dtfac/npts
	cntout = round(Int,dtout/dt)
	nsteps = round(Int,tfin/dt)
	epsilon = epsfac/npts
	sigma = sigfac/npts
	params = ParamType(dt,epsilon,sigma,0,lenevo)

	# Create the folders for saving the data and plotting figures
	datafolder = "../datafiles/run/"
	newfolder(datafolder)
	plotfolder = "../figs/"
	newfolder(plotfolder)
	paramsout = string(datafolder,"params.dat")
	# Copy the parameters file to the output folder along with some additional information.
	paramvec = readparams()
	iostream = open(paramsout, "w")
	writedlm(paramsout, [paramvec; cntout; 0.])
	close(iostream)
	# Set up the target points to measure u,v,p.
	ntar0 = 10; xmax = 2.8; ymax = 0.8
	ntar,xtar,ytar,utar,vtar,ptar = targets(ntar0,xmax,ymax)

	# Begin the erosion computation, time it, and save data.
	t0 = time()
	plotnsave(thlenvec0,datafolder,plotfolder,0)
	# Use the Runge-Kutta starter.
	thlenvec1 = RKstarter!(thlenvec0, params)
	if cntout==1; plotnsave(thlenvec1,datafolder,plotfolder,1); end
	# Enter the time loop to use the multi-step method and save the data.
	for cnt = 2:nsteps
		utar,vtar,ptar = stokes!(thlenvec1,sigma,ntar,xtar,ytar)
		advance_thetalen!(thlenvec1,thlenvec0,params)
		# Save the data.
		if mod(cnt,cntout)==0
			plotnsave(thlenvec1,datafolder,plotfolder,cnt)
		end
		# Time the computation and write it to the params file.
		t1 = time(); elapsedtime = t1-t0; updatecputime(paramsout, elapsedtime)
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
