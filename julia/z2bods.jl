# Test multiple bodies.
include("spectral.jl")
include("thetalen.jl")
include("geometries.jl")
using Winston

npts = 64
nbods = 2
dth = 2*pi/npts
theta = range(0, dth, npts)
r1 = 0.2
r2 = 0.2

ntot = npts*nbods
xc1 = 0.0
yc1 = 0.5
xc2 = 0.0
yc2 = -0.5

x1 = xc1 + r1*cos(theta)
y1 = yc1 + r1*sin(theta)
x2 = xc2 + r2*cos(theta)
y2 = yc2 + r2*sin(theta)
xx = [x1; x2]
yy = [y1; y2]
tau = zeros(Float64, ntot)
tau = stokes(npts,nbods,xx,yy)

tau1 = tau[1:npts]
tau2 = tau[npts+1:end]

p1 = plot(x1,y1, x2,y2)
xlim(-1,1); ylim(-1,1)
figure(width=400, height=400); display(p1)

p2 = plot(theta,tau1,"g.-", theta,tau2,"b.-")
figure(width=600, height=400); display(p2)