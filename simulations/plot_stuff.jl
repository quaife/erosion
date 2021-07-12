# OBJECTIVE: Plot some quantities from the post-processed erosion output data files.

using Erosion

using FileIO
using Plots; plotlyjs()

#= Backends for Plots
plotlyjs: Makes two separate displays automatically, even without the reuse=false statement.

plotly: Plots in web-browser which I find annoying.

pyplot: Only makes two separate displays with the reuse=false statement, but the default size of the 
displays is too tiny and I have to manually resize them.

gr: This gives the same type of plot as with no backend specified so it must be the default.
It will not give two separate displays, even with the reuse=false statement.
=#

function plot_stuff(infile::AbstractString)
	# Read variables
	params, thldvec, areas, resist, resist_rot = 
	load(infile, "params", "thldvec", "areas", "resist", "resist_rot")

	# Extract basic stuff.
	nlast = length(thldvec)	

	# Calculate the total area at each time step and plot it.
	areatot, tt = zeros(Float64, nlast), zeros(Float64, nlast)
	for nn = 1:nlast 
		areatot[nn] = sum(areas[nn])	# Sum over all the bodies.
		tt[nn] = thldvec[nn].tt
	end
	area_plot = plot(tt, areatot, xlabel="time", ylabel="total area", label="area")

	# Plot the resistivity.
	resist_plot = plot(tt, resist, reuse=false, xlabel="time", ylabel="resistivity", label="resist")
	
	display(area_plot)
	display(resist_plot)
end

#plot_stuff("output_data/proc_data-02-1.jld2")
plot_stuff("output_data/proc_data-20-2.jld2")
