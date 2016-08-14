# Test calling a Fortran code.
# workspace()
using Winston

nn = 64
rad = 0.05

dth = 2*pi/nn
theta = range(0, dth, nn)
xx = rad*cos(theta)
yy = rad*sin(theta)
tau = zeros(nn)

ccall((:stokessolver_, "libstokes.so"),
        Void, (Ptr{Int}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
        &nn, xx, yy, tau)

tauex = sin(theta)/rad
err = maximum(abs(tau-tauex))
println(err)
plot(theta,tau,".-k",theta,tauex,".-b")