#= OBJECTIVE: Output data for Veusz to plot figure 1, single-run data. =#
include("plot_defns.jl")

#= Note: the format of DragData object is:
DragData(pdragx, pdragy, vdragx, vdragy, umax, tau_all, atau_all, press_all)
tau and press take the form tau[1:npts, bod] =#

# Main routine to make plots in Veusz.
function main(label::AbstractString)
	# Read variables
	thldvec, areas_vec, resist, resist_rot, resist_circs, resist_circs_rot, drag_data = 
		load(proc_file(label), "thldvec", "areas_vec", "resist", "resist_rot", 
		"resist_circs", "resist_circs_rot", "drag_data")

	# Calculate time-dependent quantities.
	tt = [ thldvec[nn].tt for nn in eachindex(thldvec)]
	area = [ sum(areas_vec[nn]) for nn in eachindex(areas_vec) ]
	hdrag = [ drag_data[nn].pdragx + drag_data[nn].vdragx for nn in eachindex(drag_data)] 
	vdrag = [ drag_data[nn].pdragy + drag_data[nn].vdragy for nn in eachindex(drag_data)]
	umax = [ drag_data[nn].umax for nn in eachindex(drag_data) ]
	# Calculate more quantities.
	perm = permeability.(resist)
	anis = anisotropy.(resist, resist_rot)
	#perm_ratio.(resist, resist_circs)

	# Make text file for Veusz to plot.
	vplot([tt area resist resist_rot hdrag vdrag umax perm anis], vfolder("fig1"))
end

# Possible runs: 02-1, 20-2, 20-5, 20-8, 40-3, 40-7, 40-8
main("40-3")
