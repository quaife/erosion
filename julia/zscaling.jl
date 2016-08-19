# Determine how the shear stress scales with body size.
include("includes.jl")

npts = 32
nbods = 1
dth = 2*pi/npts
theta = range(0, dth, npts)

nrads = 15
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