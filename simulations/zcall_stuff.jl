# OBJECTIVE: Collection of commands to call various routines from the ErosionPackage.

using Erosion
using JLD2, FileIO
include("remake_data.jl")

# This is the list of big runs that were saved in the text format circa 2019.
good_runs() = [(20,2), (20,5), (20,8), (40,3), (40,7), (40,8), 
				(60,3), (60,7), (60,9), (80,4), (80,7), (80,9),
				(100,3), (100,6), (100,9)]

#--------------- ROUTINES FOR BATCH PROCESSING ---------------#
#= OBJECTIVE: Make many initial geometries in the JLD2 format.
These calls use the same parameters as the saved text input geometries with afrac = 0.6.
However, I discovered the geometries they make are different, likely due to changes in Julia versions. =#
function make_many_input_geos()
	afrac = 0.6
	runs = good_runs()
	for run in runs
		nbods, seed = run
		make_geos(nbods, afrac, seed, makeplots=false)
	end
end

# OBJECTIVE: Reformat many input geometries from the text to JLD2 files.
function remake_many_input_geos()
	runs = good_runs()
	for run in runs
		nbods, seed = run
		remake_input_geos("afrac06", nbods, seed)
	end
end

# OBJECTIVE: Reformat many output data sets from text to JLD2 files.
function remake_many_output_data()
	runs = good_runs()
	for run in runs
		nbods, seed = run
		remake_output_data(string(nbods, "-", seed))
	end
end

# OBJECTIVE: Post-process all of the remade data at once.
function post_proc_many()
	runs = good_runs()
	for run in runs
		nbods, seed = run
		file = string("output_data/raw_data-", nbods, "-", seed, ".jld2")
		post_process(file)
	end
end
#-------------------------------------------------#

#= OBJECTIVE: Make several new input geometries with seeds 10 and up
so that I do not resuse the name of already saved runs. Remember, for
a given seed, the construced geometry will be different due to changes 
in Julia. =#
function make_new_geos(nbods::Int, areafrac::Float64 = 0.6)
	for ii=10:19
		make_geos(nbods, areafrac, ii, makeplots=false)
	end
end

#--------------- CALL THE ROUTINES ---------------#

# CALL BASIC ROUTINES FOR TESTING
# make_geos(2, 0.02, 8)
# make9geos(2, 0.02)
# remake_input_geos("afrac06", 20, 2)
# remake_output_data("20-2")
# post_process("output_data/raw_data-20-2.jld2")

# CALL THE BATCH PROCESSING ROUTINES
# make_many_input_geos()
# remake_many_input_geos()
# remake_many_output_data()
# post_proc_many()

# CALL THE ROUTINE TO MAKE NEW GEOMETRIES.
#make_new_geos(60)