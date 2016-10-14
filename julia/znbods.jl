# Erode multiple bodies.
include("basic.jl")
function main()
	##### PARAMETERS #####
	# Geometry parameters.
	npts = 256
	nbods = 1
	xsm1, ysm1 = +0.0, +0.4
	xsm2, ysm2 = -0.0, -0.4
	xsm3, ysm3 = +0.4, +0.0
	xsm4, ysm4 = -0.4, -0.0	
	# For circle geometry
	rad1, rad2, rad3, rad4 = 0.2, 0.2, 0.2, 0.2
	# Evolution parameters.
	dt = 0.5e-3
	epsilon = 5*dt^(2/3)
	sigma = epsilon
	nsteps = 50
	# Misc parameters.
	beta = 0
	axlims = [1.,1.]
	# Target points
	ntargs = 11
	yend = 0.8
	ytar = collect(linspace(-yend, yend, ntargs))
	xtar = -2.8 * ones(Float64, ntargs)
	######################

	# Put the parameters in a single variable.
	params = ParamType(dt,epsilon,sigma,beta)
	# Initialize variables.
	utar,vtar,ptar = [zeros(Float64,ntargs) for ii=1:3]
	pavg = zeros(Float64,nsteps)
	# Create the initial geometries.
	thlen01 = circgeo(npts,rad1,xsm1,ysm1)
	thlen02 = circgeo(npts,rad2,xsm2,ysm2)
	thlen03 = circgeo(npts,rad3,xsm3,ysm3)
	thlen04 = circgeo(npts,rad4,xsm4,ysm4)
	# Create the vector of ThetaLenType values.
	thlenvec0 = [thlen01, thlen02, thlen03, thlen04]
	#thlenvec0 = [thlen01]
	
	# Plot the initial geometries, t=0.
	plotcurves!(thlenvec0,0; axlims=axlims)	
	# Use RK2 as a starter.
	thlenvec1 = RKstarter!(thlenvec0, params)
	# Plot the result for t=dt.
	plotcurves!(thlenvec1,1; axlims=axlims)

	# Enter the time loop.
	for cnt = 2:nsteps
		# Compute the new stress and save it.
		utar,vtar,ptar = stokes!(thlenvec1,sigma,ntargs,xtar,ytar)
		# Advance thlen forward in time using the multi-step method.
		advance_thetalen!(thlenvec1,thlenvec0,params)
		# Calculate the average pressure.
		pavg[cnt] = mean(ptar)
		# Plot the results.
		plotcurves!(thlenvec1,cnt; axlims=axlims)
		#plotpress(ytar,ptar,cnt)
	end
	#p1 = plot(pavg,".-")
	#display(p1)
end

main()
