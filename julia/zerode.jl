# Erode a body.
include("includes.jl")

npts = 32
nbods = 1
rad = 0.2

# Create the initial circular geometry.
theta,len = circgeo(npts,rad)
# Accomodate the possibility of multiple bodies for the Stokes solver.
thetas = theta
lens = [len]

HERE stokes_thl(npts,nbods,thetas,lens)

tm = 0.0
while(tm <= tfin)

end