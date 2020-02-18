# shrink_circles.jl
include("main.jl")

function shrink_circles(paramsfile::AbstractString = "params")
	nsteps = 100
	thlenden, params = startup(paramsfile)
	npts,nbods = getnvals(thlenden.thlenvec)
	# Initialize the length vector.
	lenvec0 = zeros(nbods)
	for  mm = 1:nbods
		lenvec0[mm] = thlenden.thlenvec[mm].len
	end
	# For each step, shrink the circles.
	for nn = 0:nsteps-1
		# Scale the lengths by a common factor.
		lfac = sqrt((nsteps - nn)/(nsteps))
		for mm = 1:nbods
			thlenden.thlenvec[mm].len = lfac*lenvec0[mm]
			thlenden.density = []
			thlenden.denrot = []
			thlenden.thlenvec[mm].xx = []
			thlenden.thlenvec[mm].yy = []
		end
		# Plot and save the results.
		plotnsave(nn,0.,thlenden,params)
	end
	postprocess(string("run_",paramsfile))
end
