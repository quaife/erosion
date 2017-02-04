
function run()
	include("basic.jl")

	thlenvec = readthlenfile(thlenfile)
	erosion(thlenvec)
end

run()



#tfin = 0.5
#dt = 0.05/npts
