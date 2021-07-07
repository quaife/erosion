#!/bin/bash

for i in {1,2,4,8,16,32,48,64}; do
	cat fmm_test$i.out | grep 'GMRES  time' | awk '{print $NF}' | awk '{x += $1}END{print x/3600/'"$i"'}'
done > fmm_times_gmres.out

for i in {1,2,4,8,16,32,48,64}; do
	cat fmm_test$i.out | grep cpu | grep -o [0-9][0-9.]*
done > fmm_times_wall.out
