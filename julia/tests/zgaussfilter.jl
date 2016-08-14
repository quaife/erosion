# Test gaussfilter
# workspace()
# set_zero_subnormals(true)
include("spectral.jl")
using Winston

function wrap(fx)
	return [fx; fx; fx[1]]
end

function teststuff()
	nn = 64
	sigma = 0.2
	# Calculated parameters
	aa = -pi
	bb = pi
	dx = (bb-aa)/nn
	xx = range(aa, dx, nn)
	# Smooth a function
	fx = abs(xx)
	fs = gaussfilter(fx,sigma)
	# Plot
	xext = [xx; xx+(bb-aa); 2*bb-aa]
	plot(xext, wrap(fx), "k.-", xext, wrap(fs), "b.-")
end

teststuff()