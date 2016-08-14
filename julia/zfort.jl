# Test calling a Fortran code.
# workspace()
include("spectral.jl")
include("thetalen.jl")
include("geometries.jl")
using Winston

nn = 128
rad = 1e-4

dth = 2*pi/nn
theta = range(0, dth, nn)
xx = rad*cos(theta)
yy = rad*sin(theta)
tau = stokes(nn,xx,yy)

tauex = -0.33/rad * sin(theta)
err = maximum(abs(tau-tauex))
println(err)
plot(theta,tau,".-k",theta,tauex,".-b")