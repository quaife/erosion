# MAIN GOAL: Postprocess the jld2 file to compute new quantities.
# Convention: el = 1:nbods indexes the bodies, nn indexes the timestep.

using LinearAlgebra
include("run0.jl")
include("main.jl")

#--- ROUTINES TO COMPUTE THE PRESSURE ON THE SURFACE ---#
#= Note: The purpose of these routines is to compute the pressure on the surface
of each body, rather than at a set of target points in the fluid domain. =#
# Fortran wrapper.
function compute_pressure(xx::Vector{Float64}, yy::Vector{Float64}, 
		density::Vector{Float64}, npts::Int, nbods::Int, nouter::Int, ibary::Int)
	pressure = zeros(Float64, npts*nbods)
	if nbods > 0
		ccall((:computepressure_, "libstokes.so"), Nothing,
			(Ref{Int},Ref{Int},Ref{Int},
			Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64}),
			npts, nbods, nouter, xx, yy, density, pressure)
	end
	return pressure
end
# Dispatch for ThLenDenType.
function compute_pressure(thlenden::ThLenDenType, params::ParamSet;
		fixpdrop::Bool=false, rotation::Bool=false)
	@unpack npts, nouter, ibary = params
	nbods,xv,yv,density = getnxyden(thlenden,params,fixpdrop,rotation)
	pressure = compute_pressure(xv,yv,density,npts,nbods,nouter,ibary)
	return pressure
end
#-------------------------------------------------#

#----------- ROUTINES FOR AREA, RESISTIVITY, DRAG, ETC. -----------#
# Get the normal and tangent directions.
# Convention: CCW parameterization and inward pointing normal.
function getns(theta::Vector{Float64}, rotation::Bool=false)
	# CCW tangent vector.
	if rotation == false
		sx, sy = cos.(theta), sin.(theta)
	else
		sx, sy = -sin.(theta), cos.(theta)
	end
	# Inward pointing normal vector.
	nx, ny = -sy, sx
	return sx, sy, nx, ny
end
# Compute the area of each body.
function getareas(thlenden::ThLenDenType)
	nbods = length(thlenden.thlenvec)
	areavec = zeros(Float64, nbods)
	for el = 1:nbods
		thlen = thlenden.thlenvec[el]
		xx, yy = getxy(thlen)
		sx,sy,nx,ny = getns(thlen.theta)
		npts = length(thlen.theta)
		ds = thlen.len / npts
		# Compute area in two ways to estimate error.
		areax = -dot(xx,nx)*ds
		areay = -dot(yy,ny)*ds
		area = 0.5*(areax+areay)
		reldiff = abs(areax-areay)/area
		reldiff > 1e-3 ? @warn(string("Relative error in area = ", round(reldiff,sigdigits=2))) : 0.
		areavec[el] = area
	end
	if nbods == 0 areavec = [0.0] end
	return areavec
end
# Compute the resistivity/permeability of the porous matrix.
function resistivity(thlenden::ThLenDenType, params::ParamSet, x0::Float64=2.0; rotation::Bool=false)
	# Retrieve the pressure drop and flux (assuming umax = 1)
	pdrop, qavg = getpdrop(thlenden, params, x0, rotation=rotation)
	# Calculate the total resistivity
	resist = pdrop/(2*x0*qavg)
	# If pipe flow (ibc = 0) then remove the contribution from the walls.
	if params.ibc == 0; resist = x0*(resist - 3); end
	return resist
end

# Structure for the drag output data.
mutable struct DragData
	pdragx::Float64; pdragy::Float64; vdragx::Float64; vdragy::Float64; 
	umax::Float64; tauvec::Vector{Float64}; atauvec::Vector{Float64};
end
# Compute the total drag on all of the bodies combined.
function drag(thlenden::ThLenDenType, params::ParamSet; rotation::Bool=false)
	# Get the shear stress and pressure on the set of bodies.
	# Note 1: the stress is not smoothed and absolute value is not taken.
	# Note 2: these are the values with umax = 1.
	tauvec = compute_stress(thlenden, params.nouter, params.ibary, fixpdrop=false, rotation=rotation)
	pressvec = compute_pressure(thlenden, params.nouter, params.ibary, fixpdrop=false, rotation=rotation)
	thlenvec = thlenden.thlenvec
	nbods = length(thlenden.thlenvec)
	npts = params.npts
	pdragx,pdragy,vdragx,vdragy = 0.,0.,0.,0.
	atauvec = zeros(Float64,npts*nbods)
	# Loop over all of the bodies.
	for el=1:nbods
		# Get the pressure and stress on body el.
		el1 = npts*(el-1)+1
		el2 = npts*el
		press = pressvec[el1:el2]
		tau = tauvec[el1:el2]
		# Go ahead and compute the absolute smoothed stress.
		atauvec[el1:el2] = gaussfilter(abs.(tau), params.sigma)
		# Get the tangent/normal vectors and arc length increment.
		sx, sy, nx, ny = getns(thlenvec[el].theta, rotation)
		ds = thlenvec[el].len / npts
		# Compute the pressure and viscous drag separately.
		# Note: I believe both should be plus signs due to the conventions of s and n.
		pdragx += dot(press,nx)*ds
		pdragy += dot(press,ny)*ds
		vdragx += dot(tau,sx)*ds
		vdragy += dot(tau,sy)*ds
	end
	umax = getumax(thlenden, params.nouter, params.ibary, params.fixpdrop)
	return DragData(pdragx, pdragy, vdragx, vdragy, umax, tauvec, atauvec)
end
#-------------------------------------------------#

#----------- ROUTINES FOR TARGET POINT CALCULATIONS -----------#
# Set up target points on a body-fitted grid around all bodies.
function bodyfitgrid(thlenv::Vector{ThetaLenType}, spacevec::Vector{Float64}, nptslayer::Int)
	nbods = length(thlenv)
	if nbods == 0; return [],[]; end
	npts = length(thlenv[1].theta)
	# Use nptslayer in each layer.
	ind0 = div(npts, 2*nptslayer)
	ind0 = max(ind0, 1)
	ind = ind0 : 2*ind0 : npts
	nlayers = length(spacevec)
	xtar, ytar = Array{Float64}(undef,0), Array{Float64}(undef,0)
	# Loop over the bodies.
	for el = 1:nbods
		thlen = thlenv[el]
		xx, yy = thlen.xx[ind], thlen.yy[ind]
		nx, ny = getns(thlen.theta)[3:4]
		nx, ny = nx[ind], ny[ind]
		append!(xtar, vec(xx*ones(1,nlayers) - nx*transpose(spacevec)))
		append!(ytar, vec(yy*ones(1,nlayers) - ny*transpose(spacevec)))
		# Remove the points that lie outside the computational domain.
		badind = findall( abs.(ytar) .> 0.999 )
		deleteat!(xtar, badind)
		deleteat!(ytar, badind)
	end
	return xtar, ytar
end
# regbodtargs: Set up target points on a regular and body fitted grid.
function regbodtargs(thlenv::Vector{ThetaLenType})
	# Make the regular grid.
	hh = 0.05
	xlocs = collect(-1-2*hh: hh: 1+2*hh)
	ylocs = collect(-1+0.5*hh: hh: 1-0.5*hh)	
	xreg, yreg = regulargrid(xlocs, ylocs)
	# Make the body-fitted grid.
	spacevec = 0.01*collect(1:2:3)
	nptslayer = 32
	xbod, ybod = bodyfitgrid(thlenv, spacevec, nptslayer)
	# Combine the regular and body fitted grid into a single set of points.
	targets = TargetsType([], [], [], [], [], [])
	targets.xtar = [xreg; xbod]
	targets.ytar = [yreg; ybod]
	return targets
end
#-------------------------------------------------#


#--------------------- MAIN ROUTINES ------------------#
# Read the basic parameters from the jld2 file.
function read_vars(datafile::AbstractString)
	file = jldopen(datafile, "r")
	params = read(file, "params")
	thldvec = read(file, "thldvec")
	close(file)
	return params, thldvec
end

# pp1: Postprocess the fast stuff: area and resistivity.
function pp1(datafile::AbstractString)
	println("Beginning pp1 on ", datafile)
	params, thldvec = read_vars(datafile)
	nlast = length(thldvec)
	areas, resist, resist_rot = [], [], []
	# Loop over the time values to compute areas and resistivity at each.
	for nn = 1:nlast
		println("pp1 step ", nn, " of ", nlast)
		thlenden = thldvec[nn]
		# Compute the area of each body.
		push!(areas, getareas(thlenden))		
		# Compute the resistivity and push to the output vectors.
		push!(resist, resistivity(thlenden, params))
		push!(resist_rot, resistivity(thlenden, params, rotation=true))
	end
	# Save the new data to the same jld2 file.
	jldopen(datafile, "r+") do file
		write(file, "areas", areas)
		write(file, "resist", resist)
		write(file, "resist_rot", resist_rot)
	end
	println("Finished pp1 on ", datafile, "\n")
end

# pp2: Postprocess the slower stuff: drag and stress.
function pp2(datafile::AbstractString)
	println("\n\nBeginning pp2 on ", datafile)
	params, thldvec = read_vars(datafile)
	nlast = length(thldvec)
	drag_data, drag_data_rot = [], []
	# Loop over the time values to compute the drag and stress at each.
	for nn = 1:nlast
		println("pp2 step ", nn, " of ", nlast)
		thlenden = thldvec[nn]
		# Compute the drag and stress, and push to the output vectors.
		push!(drag_data, drag(thlenden, params) )
		push!(drag_data_rot, drag(thlenden, params, rotation=true) )
	end
	# Save the new data to the same jld2 file.
	jldopen(datafile, "r+") do file
		write(file, "drag_data", drag_data)
		write(file, "drag_data_rot", drag_data_rot)
	end
	println("Finished pp2 on ", datafile, "\n")
end

# pp3: Postprocess the slowest stuff: quantities of interest at the target points.
function pp3(datafile::AbstractString)
	println("\n\nBeginning pp3 on ", datafile)
	params, thldvec = read_vars(datafile)
	nlast = length(thldvec)	
	target_data = []
	# Loop over the time values to compute the target-point data at each.
	for nn = 1:nlast
		println("pp3 step ", nn, " of ", nlast)
		thlenden = thldvec[nn]
		# Compute velocity, pressure, vorticity at a set of target points, with umax set to 1.
		targets = regbodtargs(thlenden.thlenvec)
		compute_qoi_targets!(thlenden, targets, params.nouter, params.ibary, fixpdrop=false)
		push!(target_data, targets)
	end
	# Save the new data to the same jld2 file.
	add_data(datafile, "target_data", target_data)
	println("Finished pp3 on ", datafile, "\n")
end



#postprocess: Run all postprocess routines pp1, pp2, and pp3.
function postprocess(datafile::AbstractString)
	println("\n\n%------------------------------------------------------%")
	println("Beginning postprocessing ", datafile, "\n")
	t1 = @elapsed 	pp1(datafile)
	t2 = @elapsed	pp2(datafile)
#	t3 = @elapsed	pp3(datafile)
	println("Finished postprocessing ", datafile)
	println("%------------------------------------------------------%\n\n")
end
#-------------------------------------------------#

postprocess("../output_data/data-02-1.jld2")