#!/bin/bash

# Run all models in this directory. Modify the loops below for running only a
# subset of models

processes=8
ASPECT_EXEC="../plugin/aspect"

for stokes_degree in 2 3; do # 2 3
  for discontinuous_pressure in false; do # true
    for refinement in 2 3 4 5 6 7 8; do # 2 3 4 5 6 7 8
        echo "subsection Discretization" > current.prm
        echo "  set Stokes velocity polynomial degree = $stokes_degree" >> current.prm
        echo "  set Use locally conservative discretization = $discontinuous_pressure" >> current.prm
        echo "end" >> current.prm

        echo "subsection Mesh refinement" >> current.prm
        echo "  set Initial global refinement = $refinement" >> current.prm
        echo "end" >> current.prm

        echo "set Output directory = Q${stokes_degree}_P${discontinuous_pressure}_refinement${refinement}" >> current.prm
        echo "Starting Q${stokes_degree}_P${discontinuous_pressure}_refinement${refinement}"
        cat rigid_shear.prm current.prm | mpirun -np $processes $ASPECT_EXEC --
    done
  done
done
