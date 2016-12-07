# main.jl

# erosion: The main routine to erode a group of bodies for input of Vector{ThetaLenType}.
function erosion(tfin::Float64, dt::Float64, thlenvec0::Vector{ThetaLenType}; 
	axlims::Vector{Float64} = [1.,1.])
	# Extract the basic parameters
	npts = endof(thlenvec0[1].theta)
	nbods = endof(thlenvec0)
	nsteps = round(Int,tfin/dt)
	# Calculate the smoothing parameters based on the spatial resolution.
	epsilon = 10./npts
	sigma = 15./npts
	params = ParamType(dt,epsilon,sigma,0)
	# Set up the target points to measure u,v,p.
	ntar0 = 10; xmax = 2.8; ymax = 0.8
	ntar,xtar,ytar,utar,vtar,ptar = targets(ntar0,xmax,ymax)

	# Use the Runge-Kutta starter, while plotting and saving before and after.
	plotnsave(thlenvec0,0,axlims=axlims)
	thlenvec1 = RKstarter!(thlenvec0, params)
	plotnsave(thlenvec1,1,axlims=axlims)
	# Enter the time loop.
	for cnt = 2:nsteps
		utar,vtar,ptar = stokes!(thlenvec1,sigma,ntar,xtar,ytar)
		advance_thetalen!(thlenvec1,thlenvec0,params)
		plotnsave(thlenvec1,cnt,axlims=axlims)
	end
	return
end
# targets: Set up the target points to measure velocity and pressure: u,v,p.
function targets(nn::Integer, xmax::Float64, ymax::Float64)
	# Make the grid.
	ytar = collect(linspace(-ymax,ymax,nn))
	ytar = [ytar; ytar]
	xtar = ones(Float64,nn)
	xtar = xmax*[-xtar; xtar]
	# Initialize u,v,p at target points.
	utar,vtar,ptar = [zeros(Float64,2*nn) for ii=1:3]
	return 2*nn,xtar,ytar,utar,vtar,ptar
end
# plotnsave: Calls plotcurves() and savexydata()
function plotnsave(thlenvec::Vector{ThetaLenType}, cnt::Integer;
		axlims::Vector{Float64}=[3.,1.] )
	plotfile = string("../figs/shape", string(cnt), ".pdf")
	savefile = string("../datafiles/geoout", string(cnt), ".dat")
	plotcurves(thlenvec,plotfile,axlims=axlims)
	savexydata(thlenvec,savefile)
end
