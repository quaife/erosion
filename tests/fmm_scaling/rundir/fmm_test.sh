#!/bin/bash

for i in {1,2,4,8,16,32,48,64}; do
	export LD_LIBRARY_PATH="/home/alemmon/cleanroom/cleanroom/fortran_code"
	export OMP_NUM_THREADS=$i
	/home/alemmon/cleanroom/cleanroom/rundir/julia-1.5.3/bin/julia fmm_test.jl &> fmm_test$i.out 
done
