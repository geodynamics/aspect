#!/bin/bash

# Instructions for how to use this script are provided in the README.

for stokes_degree in 2; do #3
  for discontinuous_pressure in false; do # true
    for refinement in 3 4; do
      for particles_per_direction in 4 5; do # 6 7 10 15 20
        echo "subsection Discretization" > current.prm
        echo "  set Stokes velocity polynomial degree = $stokes_degree" >> current.prm
        echo "  set Use locally conservative discretization = $discontinuous_pressure" >> current.prm
        echo "end" >> current.prm

        echo "subsection Postprocess" >> current.prm
        echo "  subsection Particles" >> current.prm
        echo "    subsection Generator" >> current.prm
        echo "      subsection Reference cell" >> current.prm
        echo "        set Number of particles per cell per direction = $particles_per_direction" >> current.prm
        echo "      end" >> current.prm
        echo "    end" >> current.prm
        echo "  end" >> current.prm
        echo "end" >> current.prm

        echo "subsection Mesh refinement" >> current.prm
        echo "  set Initial global refinement = $refinement" >> current.prm
        echo "end" >> current.prm


        echo "set Output directory = Q${stokes_degree}_P${discontinuous_pressure}_refinement${refinement}_${particles_per_direction}" >> current.prm
        echo "Starting Q${stokes_degree}_P${discontinuous_pressure}_refinement${refinement}_${particles_per_direction}"
        cat time_dependent_annulus.prm current.prm | mpirun -np 4 ./plugin/aspect --
      done
    done
  done
done
