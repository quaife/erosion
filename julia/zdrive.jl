# Driver
function driver(thlenfile::AbstractString, tfin::Float64=0.5)
	include("basic.jl")
	# Extract the geometry from a data file.
	thlenvec,npts = readthlenfile(thlenfile)
	# Optional: set dt based on npts.
	dt = 0.1/npts
	# Call the main erosion routine.
	erosion(tfin,dt,thlenvec; axlims=[1.,1.])
end

function run()
	include("basic.jl")
	thlenfile = "../datafiles/thlen.dat"
	make4circs(thlenfile,256,4)
	tfin = 0.03
	driver(thlenfile,tfin)
end
run()
