# Run a simulation!

include("erosion.jl")
#using .ErosionSimulation #: run_erosion, ParamSet

function run()
	params = ParamSet()
	println(params)

	run_erosion(params)
end

run()