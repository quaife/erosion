# main.jl: The main routines to call
using Plots
using DelimitedFiles
using Statistics
using FFTW
include("basic.jl")
##include("makegeos.jl")
include("spectral.jl")
include("thetalen.jl")
include("ioroutines.jl")


sig(var::AbstractFloat, sigdig::Int=3) = round(var,sigdigits=sigdig)


#--------------- MAIN ROUTINE ---------------#
# The main routine to erode a group of bodies.
function erosion(params::paramset)
	thlenden = get_thlenden(params)
	nn=0; nfile = 0; tt = 0.;
	plotnsave(nfile,tt,thlenden,params)
	# Enter the time loop to apply Runge-Kutta.
	while(tt < params.tfin - 0.1*params.dt && length(thlenden.thlenvec) > 0)
		# Print statements
		nn += 1
		println("\n\n\nTIME STEP ", nn)
		println("t/tfin = ", sig(tt/params.tfin, 3))
		# Advance the variables forward one timestep with RK4.
		thlenden, dt = rungekutta2(thlenden, params)
		tt += dt
		# Plot and save the data if appropriate.
		if mod(nn, params.cntout)==0
			nfile += 1
			plotnsave(nfile,tt,thlenden,params)
		end
	end
	# Plot and save one last time with zero bodies.
	nfile += 1
	plotnsave(nfile,tt,thlenden,params)
	return thlenden,params,tt
end

# Get initial thlenden.
function get_thlenden(params::paramset)
	circdata = readvec(params.infile)
	nbods = round(Int, circdata[1])
	deleteat!(circdata,1)
	thlenvec = new_thlenvec(nbods)
	for nn = 1:nbods
		rad, xc, yc = circdata[1], circdata[2], circdata[3]
		thlenvec[nn] = circ2thlen(params.npts, rad, xc, yc)
		deleteat!(circdata,1:3)
	end
	return new_thlenden(thlenvec)
end









	datafolder, plotfolder = getfoldernames(params.paramsfile)
	newfolder(datafolder)
	newfolder(plotfolder)
	# Save the input parameters file.

