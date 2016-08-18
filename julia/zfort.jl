# Test calling a Fortran code.
# workspace()
include("spectral.jl")
include("thetalen.jl")
include("geometries.jl")
using Winston

npts = 32
nbods = 1
dth = 2*pi/npts
theta = range(0, dth, npts)

nrads = 10
rav = zeros(Float64,nrads)
mtv = zeros(Float64,nrads)
rad = 0.2
for ii=1:nrads
	xx = rad*cos(theta)
	yy = rad*sin(theta)
	tau = stokes(npts,nbods,xx,yy)
	maxtau = maximum(abs(tau))
	# Save vectors
	rav[ii] = rad
	mtv[ii] = maxtau
	rad = 0.5*rad
end

# Plot to determine scaling of tau with rad.
# Compare against the law 1/(r log r).
rpl = 2./( -rav.*log(rav)  )
loglog(rav,mtv,".-",rav,rpl,".-r")
xlabel("radius"); ylabel("max tau")

#= Compare the stress against a guess of the exact solution.
tauex = -0.33/rad * sin(theta)
err = maximum(abs(tau-tauex))
println(err)
plot(theta,tau,".-k",theta,tauex,".-b")
=#