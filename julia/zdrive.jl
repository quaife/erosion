# Driver
function run()
	include("basic.jl")
	# Parameters
	thlenfile = "../datafiles/thlen.dat"
	npts = 256
	nbods = 4
	dt = 0.1/npts
	tfin = 0.03
	# Make the geometry and write to a data file.
	make4circs(thlenfile,npts,nbods)
	# Extract the geometry from a data file.
	thlenvec,npts = readthlenfile(thlenfile)
	# Call the main erosion routine.
	erosion(tfin,dt,thlenvec; axlims=[1.,1.])
end

run()
