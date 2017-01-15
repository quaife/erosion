# Driver
function driver(geoinfile::AbstractString, tfin::Float64=0.5)
	include("basic.jl")
	# Extract the geometry from a data file.
	thlenvec,npts = readgeodata(geoinfile)
	# Optional: set dt based on npts.
	dt = 0.1/npts
	# Call the main erosion routine.
	erosion(tfin,dt,thlenvec; axlims=[1.,1.])
end

function run()
	include("basic.jl")
	geoinfile = "../datafiles/thlen.dat"
<<<<<<< HEAD
	make4circs(geoinfile,256,4)
=======
# 	make4circs(geoinfile,256,4)
>>>>>>> 7aaa6d37f8148cb704db648bf89296c3ff8aa0f7
	tfin = 0.03
	driver(geoinfile,tfin)
end
run()
