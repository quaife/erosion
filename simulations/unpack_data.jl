#= OBJECTIVE: Some commands to unpack the data from the erosion simulations,
mostly to help Jake to process the tortuosity from my saved runs.
This code will unpack the post-processed data files.=#

using Erosion, FileIO, DelimitedFiles, Interpolations, Statistics
using Plots; plotlyjs()

#---------------------------------------------------------------#
# Set file and folder names.
proc_file(run::AbstractString) = "output_data/proc_data-$(run).jld2"

# Calculate the porosity versus time for a given run.
function get_porosity(run::AbstractString)
	areas_vec = load(proc_file(run), "areas_vec")
	area = [ 0.25*sum(areas_vec[nn]) for nn in eachindex(areas_vec) ]
	porosity = 1 .- area
	return porosity
end

#= Read the thldvec variable. This variable contains the theta-L 
data and the density function at each time step. The index corresponds
to the time step. =#
function do_stuff()
	run = "20-2"
	thldvec = load(proc_file(run), "thldvec")
	thlenden_circs_vec = load(proc_file(run), "thlenden_circs_vec")

	for n in eachindex(thldvec)
		thlenvec = thldvec[n].thlenvec
		density, denrot = thldvec[n].density, thldvec[n].denrot
		
		#=for bod in eachindex(thlenvec)
			theta, len = thlenvec[bod].theta, thlenvec[bod].len, 
			xsm, ysm = thlenvec[bod].xsm, thlenvec[bod].ysm
		end=#
	end
end

#do_stuff()
