# OBJECTIVE: An example file that shows how to run an erosion simulation.
# This file calls the erosion package, initates parameters, 
# runs the erosion simulation, and then post-processes the data.

using Erosion

params = ParamSet()
run_erosion(params)
post_process(params)