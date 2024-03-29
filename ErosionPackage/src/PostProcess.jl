# MAIN GOAL: Post-process the jld2 file to compute new quantities.

module PostProcess

export getns, getareas, resistivity, compute_pressure, drag, bodyfitgrid, regbodtargs, pp1, pp2, pp3, post_process

using Erosion: outfile, add_data, circ2thlen
using Erosion.ThetaLen
using Erosion.DensityStress
using Erosion.SpectralMethods
using LinearAlgebra	# Used for dot() in drag() and for the matrix calculations in bodyfitgrid().
using Parameters: @unpack
using JLD2

# Set the name of the output file with the processed data.
procfile(params::ParamSet) = string("proc_data-", params.label, ".jld2")

# Note: the libstokes path is now imported from the ThetaLen module.

#----------- ROUTINES FOR AREA AND RESISTIVITY -----------#
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
	for bod = 1:nbods
		thlen = thlenden.thlenvec[bod]
		xx, yy = getxy(thlen)
		sx,sy,nx,ny = getns(thlen.theta)
		npts = length(thlen.theta)
		ds = thlen.len / npts
		# Compute area in two ways to estimate error.
		# Subtracting the mean of x,y makes the formula more accurate.
		areax = max(-dot(xx.-thlen.xsm, nx)*ds, 0.)
		areay = max(-dot(yy.-thlen.ysm, ny)*ds, 0.)
		area = 0.5*(areax+areay)
		absdiff = abs(areax-areay)
		if absdiff > 6e-3
			@warn(string("Body ", bod, ", area error = ", round(absdiff, sigdigits=2)))
		end
		areavec[bod] = area
	end
	if nbods == 0 areavec = [0.0] end
	return areavec
end

# Compute the resistivity of the porous matrix.
# Note: resisitivity = 1/permeability.
function resistivity(thlenden::ThLenDenType, params::ParamSet, x0::Float64=2.0; rotation::Bool=false)
	# Retrieve the pressure drop and flux (assuming umax = 1).
	pdrop, qavg = getpdrop(thlenden, params, x0, rotation=rotation)
	# Calculate the total resistivity.
	# Note: I corrected this formula by removing an erroneous x0 from the denominator.
	resist = 0.5*pdrop/qavg
	# If pipe-flow BCs, correct for the wall resistance.
	params.ibc == 0 && (resist += -3*x0)
	return resist
end

#= Compute a collection of circles with same centers and areas as the eroded configuration
at a given time. =#
function get_circs(thlenden::ThLenDenType, params::ParamSet, areas::Vector{<:AbstractFloat})
	nbods = length(thlenden.thlenvec)
	thlen_circ_vec = Vector{ThetaLenType}(undef, 0)
	for bod = 1:nbods
		thlen = thlenden.thlenvec[bod]
		rad = sqrt(areas[bod]/pi)
		npts = length(thlen.theta)
		thlen_circ = circ2thlen(npts, rad, thlen.xsm, thlen.ysm)
		# ONLY SAVE THE CIRCLE IF IT HAS POSITIVE RADIUS
		if rad > eps(1.0)
			push!(thlen_circ_vec, thlen_circ)
		end
	end
	println()
	thlenden_circs = new_thlenden(thlen_circ_vec)
	# Compute the density of the new circle configuration - time intensive!
	compute_density!(thlenden_circs, params)
	compute_density!(thlenden_circs, params, rotation=true)
	return thlenden_circs
end
#----------------------------------------------------------#

#------ ROUTINES TO COMPUTE THE PRESSURE ON THE SURFACE ------#
#= Note: The purpose of these routines is to compute the pressure on the surface
of each body, rather than at a set of target points in the fluid domain. 
This pressure computation is used to compute the total drag. =#

# Compute the pressure: Fortran wrapper.
function compute_pressure(xx::Vector{Float64}, yy::Vector{Float64}, 
		density::Vector{Float64}, npts::Int, nbods::Int, nouter::Int, ibary::Int)
	pressure = zeros(Float64, npts*nbods)
	if nbods > 0
		ccall((:computepressure_, libstokes), Nothing,
			(Ref{Int},Ref{Int},Ref{Int},
			Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64}),
			npts, nbods, nouter, xx, yy, density, pressure)
	end
	return pressure
end
# Compute the pressure: Dispatch for ThLenDenType.
function compute_pressure(thlenden::ThLenDenType, params::ParamSet;
		fixpdrop::Bool=false, rotation::Bool=false)
	@unpack npts, nouter, ibary = params
	nbods,xv,yv,density = getnxyden(thlenden,params,fixpdrop,rotation)
	pressure = compute_pressure(xv,yv,density,npts,nbods,nouter,ibary)
	return reshape(pressure, npts, nbods)
end
#----------------------------------------------------------#

#------ ROUTINES TO COMPUTE THE DRAG ON THE SURFACE ------#
# Structure for the drag output data.
mutable struct DragData
	pdragx::Float64; pdragy::Float64; vdragx::Float64; vdragy::Float64; umax::Float64; 
	tau_all::Array{Float64}; atau_all::Array{Float64}; press_all::Array{Float64}
end
# Compute the total drag on all of the bodies combined.
function drag(thlenden::ThLenDenType, params::ParamSet; rotation::Bool=false)
	# Basic parameters.	
	thlenvec = thlenden.thlenvec
	nbods = length(thlenvec)
	npts = params.npts
	# Get the shear stress and pressure on the set of bodies.
	# Note 1: the stress is not smoothed and absolute value is not taken.
	# Note 2: these are the values with umax = 1.
	tau_all = compute_stress(thlenden, params, fixpdrop=false, rotation=rotation)
	press_all = compute_pressure(thlenden, params, fixpdrop=false, rotation=rotation)
	atau_all = zeros(Float64, npts, nbods)
	pdragx, pdragy, vdragx, vdragy = 0., 0., 0., 0.
	# Loop over all of the bodies.
	for bod = 1:nbods
		# Get the pressure and stress on the given body.
		press = press_all[:, bod]
		tau = tau_all[:, bod]
		# Go ahead and compute the absolute smoothed stress.
		atau_all[:, bod] = gaussfilter(abs.(tau), params.sigma)
		# Get the tangent/normal vectors and arc length increment.
		sx, sy, nx, ny = getns(thlenvec[bod].theta, rotation)
		ds = thlenvec[bod].len / npts
		# Compute the pressure and viscous drag separately.
		# Note: I believe both should be plus signs due to the conventions of s and n.
		pdragx += dot(press, nx)*ds
		pdragy += dot(press, ny)*ds
		vdragx += dot(tau, sx)*ds
		vdragy += dot(tau, sy)*ds
	end
	umax = getumax(thlenden, params, params.fixpdrop)
	return DragData(pdragx, pdragy, vdragx, vdragy, umax, tau_all, atau_all, press_all)
end
#----------------------------------------------------------#

#----------- ROUTINES FOR TARGET POINT CALCULATIONS -----------#
# Set up target points on a body-fitted grid around all bodies.
function bodyfitgrid(thlenv::Vector{ThetaLenType}, spacevec::Vector{Float64}, nptslayer::Int)
	nbods = length(thlenv)
	if nbods == 0; return [],[]; end
	npts = length(thlenv[1].theta)
	# Use nptslayer in each layer.
	idx0 = div(npts, 2*nptslayer)
	idx0 = max(idx0, 1)
	idx = idx0 : 2*idx0 : npts
	nlayers = length(spacevec)
	xtar, ytar = Array{Float64}(undef,0), Array{Float64}(undef,0)
	# Loop over the bodies.
	for bod = 1:nbods
		thlen = thlenv[bod]
		xx, yy = getxy(thlen)
		nx, ny = getns(thlen.theta)[3:4]
		xx, yy, nx, ny = xx[idx], yy[idx], nx[idx], ny[idx]
		append!(xtar, vec(xx*ones(1,nlayers) - nx*transpose(spacevec)))
		append!(ytar, vec(yy*ones(1,nlayers) - ny*transpose(spacevec)))
		# Remove the points that lie outside the computational domain.
		badidx = findall( abs.(ytar) .> 0.999 )
		deleteat!(xtar, badidx)
		deleteat!(ytar, badidx)
	end
	return xtar, ytar
end
# Set up target points on both a regular and a body fitted grid.
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
# Post-process the slowest stuff: area, resistivity, and resistivity of circular configurations.
function pp1(params::ParamSet, thldvec::Vector{ThLenDenType})
	println("Beginning pp1:")
	nlast = length(thldvec)
	areas_vec, resist, resist_rot, resist_circs, resist_circs_rot = [[] for ii = 1:5]
	thlenden_circs_vec = []
	# Loop over the time values to compute areas and resistivity at each.
	for nn = 1:nlast
		println("pp1 step ", nn, " of ", nlast)
		thlenden = thldvec[nn]	
		# Compute the area of each body.
		areas = getareas(thlenden)
		push!(areas_vec, areas)
		# Compute a configuration of circles with same centers and areas.
		thlenden_circs = get_circs(thlenden, params, areas)
		# Compute the resistivity and push to the output vectors.
		push!(resist, resistivity(thlenden, params))
		push!(resist_rot, resistivity(thlenden, params, rotation=true))
		push!(resist_circs, resistivity(thlenden_circs, params))
		push!(resist_circs_rot, resistivity(thlenden_circs, params, rotation=true))
		push!(thlenden_circs_vec, thlenden_circs)
	end
	# Save the new data to the same jld2 file.
	jldopen(procfile(params), "r+") do file
		write(file, "areas_vec", areas_vec)
		write(file, "resist", resist)
		write(file, "resist_rot", resist_rot)
		write(file, "resist_circs", resist_circs)
		write(file, "resist_circs_rot", resist_circs_rot)
		write(file, "thlenden_circs_vec", thlenden_circs_vec)
	end
	println("Finished pp1.\n")
end

# Post-process fast stuff: drag and stress.
function pp2(params::ParamSet, thldvec::Vector{ThLenDenType})
	println("\n\nBeginning pp2:")
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
	jldopen(procfile(params), "r+") do file
		write(file, "drag_data", drag_data)
		write(file, "drag_data_rot", drag_data_rot)
	end
	println("Finished pp2.\n")
end

# Post-process fast stuff: quantities of interest at the target points.
function pp3(params::ParamSet, thldvec::Vector{ThLenDenType})
	println("\n\nBeginning pp3:")
	nlast = length(thldvec)	
	target_data = []
	# Loop over the time values to compute the target-point data at each.
	for nn = 1:nlast
		println("pp3 step ", nn, " of ", nlast)
		thlenden = thldvec[nn]
		# Compute velocity, pressure, vorticity at a set of target points, with umax set to 1.
		targets = regbodtargs(thlenden.thlenvec)
		compute_qoi_targets!(thlenden, targets, params, fixpdrop=false)
		push!(target_data, targets)
	end
	# Save the new data to the same jld2 file.
	add_data(procfile(params), "target_data", target_data)
	println("Finished pp3.\n")
end

# Run all post-process routines pp1, pp2, and pp3.
function post_process(infile::AbstractString)
	println("\n\n%------------------------------------------------------%")
	println("Beginning post-processing ", infile, "\n")
	
	# Read the variables from the raw data file.
	params, thldvec, cpu_hours = load(infile, "params", "thldvec", "cpu_hours")
	println(params)
	# Initialize the processed data file by saving the basic variables there.
	jldsave(procfile(params); params, thldvec, cpu_hours)

	# Call the post-processing subroutines.
	t1 = @elapsed 	pp1(params, thldvec)
	t2 = @elapsed	pp2(params, thldvec)
	t3 = @elapsed	pp3(params, thldvec)

	# Simplify the cpu times and print some statements.
	tmins(tt::Float64) = round(tt/60, sigdigits=2)
	println("Finished post-processing ", infile)
	println("Simulation time (hours): ", cpu_hours)
	println("Post-processing times (mins): p1 = ", tmins(t1), "; p2 = ", tmins(t2), "; p3 = ", tmins(t3))
	println("%------------------------------------------------------%\n\n")
end
# Dispatch for ParamSet so that outfile does not need to be exported.
# Use immediately after running a simulation before output files are moved.
post_process(params::ParamSet) = post_process(outfile(params))

end