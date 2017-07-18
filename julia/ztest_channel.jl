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
	r1 = resistivity(thlenden,nouter,1.2)
	r2 = resistivity(thlenden,nouter,1.5)
	r3 = resistivity(thlenden,nouter,1.8)
end

test()