# Driver
function driver(thlenfile::AbstractString, tfin::Float64=0.5)
	include("basic.jl")
	# Extract the geometry from a data file.
	thlenvec,npts = readthlenfile(thlenfile)
	# Optional: set dt based on npts.
	dt = 0.05/npts
	# Call the main erosion routine.
	erosion(tfin,dt,thlenvec;lenevo=0,axlims=[1.,1.])
end

function run()
	include("basic.jl")
	thlenfile = "../datafiles/thlen.dat"
	npts = 1024
	rad = 0.2
	make1circ(thlenfile,npts,rad)
	tfin = 0.5
	driver(thlenfile,tfin)
end
run()
