# MAIN GOAL: Postprocess the jld2 file to compute new quantities.
# Convention: nn indexes the timestep, bod = 1:nbods indexes the bodies.

module Postprocessing

	export getns, getareas, resistivity, compute_pressure, drag, bodyfitgrid, regbodtargs, pp1, pp2, pp3, postprocess

	import Erosion.add_data
	using Erosion.ThetaLen
	using Erosion.DensityStress
	using Erosion.SpectralMethods
	using LinearAlgebra	# Used for dot() in drag() and for the matrix calculations in bodyfitgrid().
	using Parameters: @unpack
	using JLD2

	# Set the name of the output file with the processed data.
	procfile(params::ParamSet) = string("../proc_data-", params.label, ".jld2")
	
	# path to libstokes.so
	const libstokes=abspath(joinpath(@__DIR__, "../..", "fortran_code/libstokes.so"))

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
			areax = -dot(xx,nx)*ds
			areay = -dot(yy,ny)*ds
			area = 0.5*(areax+areay)
			reldiff = abs(areax-areay)/area
			reldiff > 1e-3 ? @warn(string("Relative error in area = ", round(reldiff,sigdigits=2))) : 0.
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
		resist = pdrop/(2*x0*qavg)
		# If pipe flow (ibc = 0) then remove the contribution from the walls.
		if params.ibc == 0; resist = x0*(resist - 3); end
		return resist
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
	# Postprocess the fast stuff: area and resistivity.
	function pp1(params::ParamSet, thldvec::Vector{ThLenDenType})
		println("Beginning pp1:")
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
		jldopen(procfile(params), "r+") do file
			write(file, "areas", areas)
			write(file, "resist", resist)
			write(file, "resist_rot", resist_rot)
		end
		println("Finished pp1.\n")
	end

	# Postprocess the slower stuff: drag and stress.
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

	# Postprocess the slowest stuff: quantities of interest at the target points.
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

	# Run all postprocess routines pp1, pp2, and pp3.
	function postprocess(infile::AbstractString)
		println("\n\n%------------------------------------------------------%")
		println("Beginning postprocessing ", infile, "\n")
		
		# Read the variables from the raw data file.
		file = jldopen(infile, "r")
		params = read(file, "params")
		thldvec = read(file, "thldvec")
		cpu_hours = read(file, "cpu_hours")
		close(file)
		# Initialize the processed data file by saving the basic variables there.
		jldsave(procfile(params); params, thldvec, cpu_hours)

		# Call the postprocessing subroutines.
		t1 = @elapsed 	pp1(params, thldvec)
		t2 = @elapsed	pp2(params, thldvec)
		t3 = @elapsed	pp3(params, thldvec)

		# Simplify the cpu times and print some statements.
		tmins(tt::Float64) = round(tt/60, sigdigits=2)
		println("Finished postprocessing ", infile)
		println("Simulation time (hours): ", cpu_hours)
		println("Post-processing times (mins): p1 = ", tmins(t1), "; p2 = ", tmins(t2), "; p3 = ", tmins(t3))
		println("%------------------------------------------------------%\n\n")
	end

end