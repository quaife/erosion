# shrink_circles.jl
include("main.jl")

function shrink_circles(paramsfile::AbstractString = "params")
	thlenden, params = startup(paramsfile)
	nsteps = 100
	for nn = 0:nsteps
		# Scale the lengths by a common factor.
		lfac = sqrt((nsteps - nn)/nsteps)
		npts,nbods = getnvals(thlenden.thlenvec)
		for mm = 1:nbods
			thlenden.thlenvec[mm].len *= lfac
		end
		getstress!(thlenden, params)
		# Plot and save the data.
		plotnsave(nn,0.,thlenden,params)
	end
	postprocess(string("run_",paramsfile))
end

