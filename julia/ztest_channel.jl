# Test
include("main.jl")

function test()
	nouter = 1024
	# ParamType(dt,eps,sigma,nouter,iffm,fixarea)
	params = ParamType(0., 0., 0., nouter, 1, 0)
	thlenden = ThLenDenType([],[],[])
	compute_density!(thlenden,params)

	#--------------------------------------#
	# Compute the permeability
	println()
	rt1,rm1 = resistivity(thlenden,nouter,1.5)
	rt2,rm3 = resistivity(thlenden,nouter,2.0)
	rt2,rm3 = resistivity(thlenden,nouter,2.5)
end

test()