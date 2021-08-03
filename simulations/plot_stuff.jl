# OBJECTIVE: Plot some quantities from the post-processed erosion output data files.

using Erosion, FileIO
using Plots; plotlyjs()

#= Backends for Plots
plotlyjs: Makes two separate displays automatically, even without the reuse=false statement.

plotly: Plots in web-browser which I find annoying.

pyplot: Only makes two separate displays with the reuse=false statement, but the default size of the 
displays is too tiny and I have to manually resize them.

gr: This gives the same type of plot as with no backend specified so it must be the default.
It will not give two separate displays, even with the reuse=false statement.
=#

proc_file(label::AbstractString) = string("output_data/proc_data-", label, ".jld2")
plot_folder() = "zPlots/"

function plot_stuff(label::AbstractString)
	# Read variables
	params, thldvec, areas_vec, resist, resist_rot, resist_circs = 
	load(proc_file(label), "params", "thldvec", "areas_vec", "resist", "resist_rot", "resist_circs")

	# Calculate the total area at each time step.
	nlast = length(thldvec)	
	areatot, tt = zeros(Float64, nlast), zeros(Float64, nlast)
	for nn = 1:nlast 
		areatot[nn] = sum(areas_vec[nn])	# Sum over all the bodies.
		tt[nn] = thldvec[nn].tt
	end

	# Calculate the permeability.
	function perm(resist)
		abs(resist) > 10 ?  value = 1/resist : value = 0.1
		return value
	end
	# Calculate the permeability ratios.
	prat = resist_circs ./ resist
	anis = resist_rot ./ resist
	
	# Area plot
	area_plot = plot(tt, areatot, label="area", xlabel="time", ylabel="total area")
	
	# Resistivity Plot
	resist_plot = plot(tt, resist, label="erosion", xlabel="time", ylabel="resistivity", yaxis=:log, reuse=false)
	plot!(resist_plot, tt, resist_circs, label="circles")
	plot!(resist_plot, tt, resist_rot, label="rotated")

	# Permeability Plot.
	perm_plot = plot(tt, perm.(resist), label="erosion", xlabel="time", ylabel="permeability", yaxis=:log, reuse=false)
	plot!(resist_plot, tt, resist_circs, label="circles")
	plot!(perm_plot, tt, perm.(resist_circs), label="circles")
	plot!(perm_plot, tt, perm.(resist_rot), label="rotated")
	
	# Permeability Ratio.
	prat_plot = plot(tt, prat, label="eroded to circles", xlabel="time", ylabel="permeability ratio", reuse=false)
	plot!(prat_plot, tt, anis, label="anistropy")

	# Output plots.
	#display(area_plot); display(resist_plot)
	Erosion.new_folder(plot_folder())
	savefig(area_plot, string(plot_folder(),"area.pdf") )
	savefig(resist_plot, string(plot_folder(),"resist.pdf"))
	savefig(perm_plot, string(plot_folder(),"perm.pdf"))
	savefig(prat_plot, string(plot_folder(),"prat.pdf"))
end

# Possible runs: 02-1, 20-2, 20-5, 20-8, 40-3, 40-7, 40-8
plot_stuff("40-3")
