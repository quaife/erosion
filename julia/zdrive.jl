# Driver
function driver(geoinfile::AbstractString, tfin::Float64=1.)
	include("basic.jl")
	# Extract the geometry from a data file.
	thlenvec,npts = readgeodata(geoinfile)
	# Optional: set dt based on npts.
	dt = 0.2/npts

	tfin = 7*dt
	
	# Call the main erosion routine.
	erosion(tfin,dt,thlenvec; axlims=[1.,1.])
end

function run()
	include("basic.jl")
	geoinfile = "../datafiles/thlen.dat"
	make4circs(geoinfile,64,4)
	driver(geoinfile)
end
run()