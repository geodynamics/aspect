#!/bin/bash

# number of processors to run with
nproc=8

# Benchmark cases to run
for case in 1a 1b 1c 2a 2b; do
  # Refinement levels to run, note that 6 and 7 
  # take considerable time to reach steady-state
  for refinement in 3 4 5 6 7; do
		echo "subsection Mesh refinement" > current.prm
		echo "  set Initial global refinement = ${refinement}" >> current.prm
		echo "end" >> current.prm
		echo "" >> current.prm
		echo "set Output directory=output-case${case}_ref${refinement}" >> current.prm

    echo "Starting Case $case Refinement $refinement"
    cat base_case${case}.prm current.prm | mpirun -np $nproc ./plugin/aspect -- >log.txt 2>log.txt
  done
done

rm current.prm
