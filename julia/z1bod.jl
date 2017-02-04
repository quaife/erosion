# Driver
function run()
	include("basic.jl")
	# Parameters
	npts = 256
	rad = 0.2
	tfin = 0.5
	dt = 0.05/npts
	thlenfile = "../datafiles/thlen.dat"
	# Make the geometry and write to a data file.
	make1circ(thlenfile,npts,rad)
	# Extract the geometry from the data file.
	thlenvec,npts = readthlenfile(thlenfile)
	# Call the main erosion routine.
	erosion(tfin,dt,thlenvec;lenevo=0)
end

run()
