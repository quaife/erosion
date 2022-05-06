#= OBJECTIVE: Define functions to be used for plotting routines. =#

using Erosion, FileIO, DelimitedFiles, Interpolations, Statistics
using Plots; plotlyjs()

# Set file and folder names, t-grid.
proc_file(label::AbstractString) = "output_data/proc_data-$(label).jld2"
plot_folder() = "zPlots/"	# Plot folder
vfolder(name::AbstractString) = "veusz_data/$(name).txt"
tref() = 0:0.01:0.97		# Reference t-grid

# Write a text file for Veusz to read for plots.
function vplot(data::Array, file::AbstractString)
	iostream = open(file, "w")
	writedlm(iostream, Float32.(data))
	close(iostream)
end

# Calculate the permeability from resistivity.
function permeability(resist::AbstractFloat)
	rthresh = 10.0
	resist > rthresh ? value = 1/resist : value = 1/rthresh
	return value
end
# Calculate the anistropy from resistivities.
anisotropy(resist::AbstractFloat, resist_rot::AbstractFloat) = resist_rot/resist
# Calculate the ratio of permeabilities, comparing the eroded geometry to the circular one.
perm_ratio(resist::AbstractFloat, resist_circs::AbstractFloat) = resist_circs/resist

# Interpolate a variable (e.g. resistivity) on the reference t-grid.
function interp(variable::Vector{<:AbstractFloat})
	torig = range(0, 1, length=length(variable))
	interpolant = CubicSplineInterpolation(torig, variable)
	return interpolant(tref())
end