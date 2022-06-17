#= OBJECTIVE: Output data for Veusz to plot figures. 
The first routine outputs data from a single run.
The second routine outputs data from an ensemble of runs.=#

using Erosion, FileIO, DelimitedFiles, Interpolations, Statistics
using Plots; plotlyjs()
using Infiltrator 

#---------------------------------------------------------------#
# Set file and folder names.
proc_file(run::AbstractString) = "output_data/proc_data-$(run).jld2"
tort_file(run::AbstractString) = "proc_data-tort/proc_data-$(run)_tortuosity_partial.jld2"
plot_folder() = "zPlots/"	# Plot folder
vfolder(name::AbstractString) = "veusz_data/$(name).txt"

# Write a text file for Veusz to read for plots.
function vdata(data::Array, desc::AbstractString, file::AbstractString)
	iostream = open(file, "w")
	writedlm(iostream, ["descriptor $desc"])
	writedlm(iostream, Float32.(data))
	close(iostream)
end

# Calculate the permeability from resistivity.
function permeability(resist::AbstractFloat)
	rthresh = 1.0
	resist > rthresh ? perm = 1/resist : perm = 1/rthresh
	return perm
end

# Calculate the porosity versus time for a given run.
function get_porosity(run::AbstractString)
	areas_vec = load(proc_file(run), "areas_vec")
	area = [ 0.25*sum(areas_vec[nn]) for nn in eachindex(areas_vec) ]
	porosity = 1 .- area
	return porosity
end

# Append each tortuosity data vector with ones after bodies vanish.
function append_torts!(porosity, vars...)
	for var in vars
		append!(var, ones(length(porosity)-length(var)))
	end
end

# Calculate the tortuosity anisotropy.
function tort_anisotropy(xtort, ytort)
	return (ytort .-1) ./ (xtort .-1)
end
#---------------------------------------------------------------#


#---------------------------------------------------------------#
# SINGLE RUN: Make plots in Veusz for a single simulation run.

#= Note: the format of DragData object is:
DragData(pdragx, pdragy, vdragx, vdragy, umax, tau_all, atau_all, press_all)
tau and press take the form tau[1:npts, bod] =#

function vplot_single(run::AbstractString)
	# Read variables
	thldvec, resist, resist_rot, resist_circs, resist_circs_rot, drag_data = 
		load(proc_file(run), "thldvec", "resist", "resist_rot", 
		"resist_circs", "resist_circs_rot", "drag_data")

	# Read tortuosity variables file.
	xtort, ytort, xcirctort, ycirctort = load(tort_file(run), 
		"xtortuosity", "ytortuosity", "xcirctortuosity", "ycirctortuosity")
	# tort_file() also contains the variables porosities and circporosities.

	# Calculate time-dependent quantities.
	tt = [ thldvec[nn].tt for nn in eachindex(thldvec)]
	porosity = get_porosity(run)
	area = 1 .- porosity
	hdrag = [ drag_data[nn].pdragx + drag_data[nn].vdragx for nn in eachindex(drag_data)] 
	vdrag = [ drag_data[nn].pdragy + drag_data[nn].vdragy for nn in eachindex(drag_data)]
	umax = [ drag_data[nn].umax for nn in eachindex(drag_data) ]
	
	# Modify values for which bodies are gone.
	idx = findall(area .< 1e-5)
	resist[idx] .= resist_rot[idx] .= 1e-5
	resist_circs[idx] .= resist_circs_rot[idx] .= 1e-5
	umax[idx] .= NaN
	append_torts!(porosity, xtort, ytort, xcirctort, ycirctort)

	# Calculate the anisotropy of permeability.
	anis = resist_rot ./ resist
	anis_config = resist_circs_rot ./ resist_circs
	anis_shape = anis ./ anis_config

	# Calculate the anisotropy of tortuosity.
	anis_tort = tort_anisotropy(xtort, ytort)
	anis_tort_config = tort_anisotropy(xcirctort, ycirctort)
	anis_tort_shape = anis_tort ./ anis_tort_config

	# Make text file for Veusz to plot.
	vdata([tt area porosity  umax], "time area porosity umax", vfolder("basic_vars"))
	vdata(	[hdrag vdrag resist resist_rot resist_circs resist_circs_rot], 
			"hdrag vdrag resist resist_rot resist_circs resist_circs_rot", vfolder("resist_vars"))
	vdata([anis anis_config anis_shape], "anis anis_config anis_shape", vfolder("anis_vars"))
	
	# Make Veusz text file for tortuosity.
	vdata([xtort ytort xcirctort ycirctort anis_tort anis_tort_config anis_tort_shape], 
		  "xtort ytort xcirctort ycirctort anis_tort anis_tort_config anis_tort_shape", vfolder("tort_vars"))
end

# Possible runs: 20:2,5,8; 40:3,7,8; 60:3,7,9; 80:4,7,9; 100:3,6,9
# 100-9 looks quite good showing the configuration anisotropy near the end.
# 80-4 is also great for configurational anisotropy.
# 60-3 is a nice visual run for showing channelization.
vplot_single("100-9")
#---------------------------------------------------------------#



#---------------------------------------------------------------#
# STATISTICS: Make plots in Veusz for an ensemble of runs.

# For a given number of bodies, return the labels for the saved runs.
# Possible runs: 20:2,5,8; 40:3,7,8; 60:3,7,9; 80:4,7,9; 100:3,6
function get_runs(nbods::Integer)
	if nbods == 20
		return ["20-2", "20-5", "20-8"]
	elseif nbods == 40
		return ["40-3", "40-7", "40-8"]
	elseif nbods == 60
		return ["60-3", "60-7", "60-9"]
	elseif nbods == 80
		return ["80-4", "80-7", "80-9"]
	elseif nbods == 100
		return ["100-3", "100-6", "100-9"]
	else error("nbods must be 20, 40, 60, 80, or 100")
	end
end

# Reference porosity-grid
por_grid() = range(0.401, stop=0.999, length=60)

# Interpolate a variable (e.g. resistivity) on the reference porosity-grid.
function interp(porosity, var)
	interpolant = LinearInterpolation(Float64.(porosity[1:end]), Float64.(var[1:end]))
	return interpolant(por_grid())
end

#= Read a variable from an ensemble of runs, interpolate the variable 
on the reference porosity grid, then return the ensemble. 
varname can be "resist", "resist_rot", "resist_circ", ... =#
function read_var_ensemble(file, varname, runs)
	var_ensemble = Array{Float64}(undef, length(por_grid()), 0)
	for run in runs
		var = load(file(run), varname)
		porosity = get_porosity(run)
		# Fill the remainder of var with 1 to match size with porosity.
		append!(var, ones(length(porosity)-length(var)))
		var_interp = interp(porosity, var)
		var_ensemble = [var_ensemble var_interp]
	end
	return var_ensemble
end

# Compute the mean and standard deviation of some variables and save in an array.
function mean_std(vars...)
	data = Array{Float64}(undef, length(por_grid()), 0)
	for var in vars
		vmean = mean(var, dims = 2)
		vstd = std(var, dims = 2)
		data = [data vmean vstd]
	end
	return data
end

# Enforce a lower bound on the resistance
function resist_thresh(vars...)
	for var in vars
		var = max.(var,1e-5)
	end
end

# Plot the statistics
function vplot_stats()
	nbod_list = [20 40 60 80 100] #[20, 40, 60, 80, 100]
	# For each number of bodies, collect the statistics.
	for nbod in nbod_list
		runs = get_runs(nbod)
		# Read the resistance data from an ensemble of runs.
		resist = read_var_ensemble(proc_file, "resist", runs)
		resist_rot = read_var_ensemble(proc_file, "resist_rot", runs)
		resist_circs = read_var_ensemble(proc_file, "resist_circs", runs)
		resist_circs_rot = read_var_ensemble(proc_file, "resist_circs_rot", runs)
		resist_thresh(resist, resist_rot, resist_circs, resist_circs_rot)

		# Read the tortuosity data from an ensemble of runs.
		xtort = read_var_ensemble(tort_file, "xtortuosity", runs)
		ytort = read_var_ensemble(tort_file, "ytortuosity", runs)
		xcirctort = read_var_ensemble(tort_file, "xcirctortuosity", runs)
		ycirctort = read_var_ensemble(tort_file, "ycirctortuosity", runs)

		# Calculate the permeability for the ensemble.
		perm = permeability.(resist)
		perm_rot = permeability.(resist_rot)
		perm_circs = permeability.(resist_circs)
		perm_circs_rot = permeability.(resist_circs_rot)
		
		# Calculate the permeability anistropy for the ensemble.
		anis = resist_rot ./ resist
		anis_config = resist_circs_rot ./ resist_circs
		anis_shape = anis ./ anis_config
		
		# Calculate the tortuosity anistropy for the ensemble.
		anis_tort = tort_anisotropy(xtort, ytort)
		anis_tort_config = tort_anisotropy(xcirctort, ycirctort)
		anis_tort_shape = anis_tort ./ anis_tort_config

		# Compute the mean and standard deviation, then save the data.
		#data = mean_std(perm, perm_rot, perm_circs, perm_circs_rot, anis, anis_config, anis_shape)
		data = [por_grid() mean_std(perm, perm_rot, xtort, ytort, 
			anis, anis_config, anis_shape, anis_tort, anis_tort_config, anis_tort_shape)]
		desc = string("por_grid perm$nbod,+- perm_rot$nbod,+- xtort$nbod,+- ytort$nbod,+-", 
			" anis$nbod,+- anis_config$nbod,+- anis_shape$nbod,+-",
			" anis_tort$nbod,+- anis_tort_config$nbod,+- anis_tort_shape$nbod,+-")
		vdata(data, desc, vfolder("stats$nbod"))
	end
end

vplot_stats()
