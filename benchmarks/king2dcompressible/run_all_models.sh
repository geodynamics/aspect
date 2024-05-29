#!/bin/bash

# number of processors to run with
nproc=8

for formulation in ala tala; do
  for di in 0.25 0.5 1.0; do
    for rayleigh in 1e4 1e5; do
      # Refinement levels to run, note that 6 and 7 
      # take considerable time to reach steady-state
      for refinement in 3 4 5 6 7; do
        echo "subsection Adiabatic conditions model" > current.prm
        echo "  subsection Function" >> current.prm
        echo "    set Function constants  = Di=$di, gamma=1.0" >> current.prm
        echo "  end" >> current.prm
        echo "end" >> current.prm

        echo "subsection Material model" >> current.prm
        echo "  subsection King model" >> current.prm
        echo "    set Di = $di" >> current.prm
        echo "    set Ra = $rayleigh" >> current.prm
        if [ $formulation == "tala" ]; then 
          echo "    set Use TALA = true" >> current.prm
        fi
        echo "  end" >> current.prm
        echo "end" >> current.prm

        echo "subsection Mesh refinement" >> current.prm
        echo "  set Initial global refinement = $refinement" >> current.prm
        echo "end" >> current.prm

        echo "set Output directory = Di${di}_Ra${rayleigh}_refinement${refinement}_${formulation}" >> current.prm
        echo "Starting Di${di}_Ra${rayleigh}_refinement${refinement}_${formulation}"
        cat ala.prm current.prm | mpirun -np $nproc ./plugin/aspect -- >log.txt 2>log.txt
      done
    done
  done
done

rm current.prm
