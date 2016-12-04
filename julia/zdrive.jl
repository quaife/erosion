# Driver
include("basic.jl")

function driver(tfin::Float64=1.)
	# Extract the geometry from a data file.
	thlenvec,npts = readgeodata("thlen.dat")
	# Optional: set dt based on npts.
	dt = 0.2/npts

	tfin = 7*dt
	
	# Call the main erosion routine.
	erosion(tfin,dt,thlenvec; axlims=[1.,1.])
end

make4circs("thlen.dat",64,4)
driver()