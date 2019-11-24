#!/bin/bash
for i in $(seq 1 1 100)
do
   	mpirun -np 4 ./check_memory>>ftiemposc16.txt
done
