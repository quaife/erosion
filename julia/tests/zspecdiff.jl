# Test specdiff() and specint()
# workspace()
include("spectral.jl")

function teststuff()
	nn = 4
	aa = 0.0
	bb = 2*pi
	dx = (bb-aa)/nn
	xx = range(aa, dx, nn)

	# This function is resolved exactly for nn>=4.
	a0,a1,a2,b1 = randn(4)
	fx = 0.5*a0 + a1*cos(xx) + a2*cos(2*xx) + b1*sin(xx)
	# Exact derivative.
	dfex = -a1*sin(xx) -2*a2*sin(2*xx) + b1*cos(xx)
	# Exact antiderivative, forgetting about a0 term.
	Fxex = a1*sin(xx) + 0.5*a2*sin(2*xx) - b1*cos(xx)

	# Numerical derivative using specdiff and adjusting for the interval size.
	df = specdiff(fx, bb-aa)
	# Numerical antiderivative using specint and adjusting for the interval size.
	Fx = specint(fx, bb-aa)

	# Print the errors.
	errdf = maxabs(dfex-df)
	errFx = maxabs(Fxex-Fx)
	println("specdiff error = ", errdf)
	println("specint error = ", errFx)
end

teststuff()