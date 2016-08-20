# Erode a body.
include("includes.jl")

function RKstarter()
    
end

npts = 32
nbods = 1
rad = 0.2

# Create the initial circular geometry.
theta,len = circgeo(npts,rad)
# Accomodate the possibility of multiple bodies for the Stokes solver.
thetas = theta
lens = [len]
# Make the center at the origin.
xcs = [0.0]
ycs = [0.0]

# Calculate the stress
atau0 = stokes_thl(npts,nbods,thetas,lens,xcs,ycs)


tm = 0.0
while(tm <= tfin)

end