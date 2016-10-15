# misc.jl

#################### Converting between x,y and theta,len ####################
#= get_thlen: Given the x and y coordinates, calculate theta and len. =#
function get_thlen(xx::Vector{Float64}, yy::Vector{Float64})
	# The partial derivatives dx/dalpha and dy/dalpha.
	dxda = specdiff(xx)
	dyda = specdiff(yy)
	# Calculate theta and len.
	theta = atan2(dyda,dxda)


	# HERE
	theta[theta .< theta[1]] += 2*pi



	lvec = sqrt(dxda.^2+dyda.^2)
	len = mean(lvec)
	# Check that the input is nearly equally spaced in arclength.
	relerr = maxabs(lvec-len)/len
	if relerr>1.e-2; error("The coordinates are not equally spaced"); return; end 
	return theta, len
end
#= getxy: Given theta and len, reconstruct the x and y coordinates of a body.
xsm and ysm are the boundary-averaged values.
While we're at it, we can also calculate the normal direcations. =#
function getxy(theta::Vector{Float64}, len::Float64, xsm::Float64, ysm::Float64)
	# The partial derivatives dx/dalpha and dy/dalpha.
	dx = len * (cos(theta) - mean(cos(theta)))
	dy = len * (sin(theta) - mean(sin(theta)))
	# Integrate to get the x,y coordinates; result will have mean zero.
	xx = specint(dx)
	yy = specint(dy)
	# Move to have the correct average values.
	xx += xsm
	yy += ysm
	return xx,yy
end
# getxy!: Dispatch for input of type ThetaLenType. Only computes if they are not loaded.
function getxy!(thlen::ThetaLenType)
	# Only compute xx and yy if they are not already loaded in thlen.
	if thlen.xx==[] || thlen.yy==[]
		thlen.xx, thlen.yy = getxy(thlen.theta, thlen.len, thlen.xsm, thlen.ysm)
	end
	return
end
# getnormals: Calculate the surface normals.
function getnormals(theta::Vector{Float64})
	nx = -sin(theta)
	ny = cos(theta)
	return nx, ny 
end
#################### Plotting routines ####################
# plotcurve: Plot multiple curves from the theta-len values.
function plotcurves!(thlenvec::Vector{ThetaLenType}, cnt::Integer,
		xtar::Vector{Float64}=evec(), ytar::Vector{Float64}=evec(), 
		utar::Vector{Float64}=evec(), vtar::Vector{Float64}=evec(); 
		axlims::Vector{Float64}=[3.0,1.0] )
	# Make figure of given height and preserve the aspect ratio.
	height = 400
	width = axlims[1]/axlims[2]*height
	# Make plot with given axis limits.
	p1 = plot()
	xlim(-axlims[1],axlims[1]); ylim(-axlims[2],axlims[2])
	for ii = 1:endof(thlenvec)
		thlen = thlenvec[ii]
		# Compute the xy coordinates if they are not already loaded in thlen.
		getxy!(thlen)
		xx, yy = thlen.xx, thlen.yy
		# Plot the curves.
		p1 = oplot(xx,yy,"-")
	end
	# Save the figures in a folder.
	figname = string("../figs/shape",string(cnt),".pdf")
	savefig(p1, figname, width=width, height=height)
	return
end
# plotpress: Plot the pressure distribution
function plotpress(ytar::Vector{Float64}, ptar::Vector{Float64}, 
		cnt::Integer, pmax::Float64=15.0)
	p1 = plot(ytar, ptar, ".-")
	xlim(-1., 1.); ylim(0., pmax)
	figname = string("../figs/press",string(cnt),".pdf")	
	savefig(p1, figname, width=400, height=400)
end
#################### Geometry routines ####################
# getalpha: Calculate the parameterization variable, alpha = s/L.
function getalpha(npts::Integer)
	# alpha = s/L is the parameterization variable.
	dalpha = 1.0/npts
	# Use an offset grid.
	alpha = collect(range(0.5*dalpha, dalpha, npts))
	return alpha
end
# circgeo: Creates a circle.
function circgeo(npts::Integer, rad::Float64, xsm::Float64=0.0, ysm::Float64=0.0)
	# Create a new ThetaLenType variable.
	thlen = new_thlen()
	# alpha = s/L is the parameterization variable.
	alpha = getalpha(npts)
	# theta is the tangent angle.
	thlen.theta = 0.5*pi + 2*pi*alpha
	# len is the total arclength.
	thlen.len = 2*pi*rad
	# Save xsm and ysm too.
	thlen.xsm = xsm; thlen.ysm = ysm
	return thlen
end
# polygongeo: Creates a polygon with number of sides, nsides.
function polygongeo(npts::Integer, nsides::Integer, 
		sigma::Float64 = 0.1, sdlen::Float64=0.5, xsm::Float64=0.0, ysm::Float64=0.0)
	# Create a new ThetaLenType variable.
	thlen = new_thlen()
	# alpha = s/L is the parameterization variable.
	alpha = getalpha(npts)
	# theta is the tangent angle.
	theta = zeros(npts)
	for nn=1:nsides
		a0 = (nn-1)/nsides
		a1 = nn/nsides
		intvl = ((a0 .<= alpha) & (alpha .<= a1))
		theta[intvl] = 0.5*pi + 2*pi*(nn-1)/nsides
	end
	# Smooth theta
	theta = gaussfilter(theta - 2*pi*alpha, sigma) + 2*pi*alpha
	# Save theta in the ThetaLenType object.
	thlen.theta = theta
	# len is the total arclength.
	thlen.len = nsides*sdlen
	# Save xsm and ysm too.
	thlen.xsm = xsm; thlen.ysm = ysm
	return thlen
end
##################################################
