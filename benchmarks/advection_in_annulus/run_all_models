#!/bin/bash

# This script is used to run several models consecutively,
# using all possible parameter combinations.

for temperature_degree in 1 2 3; do
  for refinement in 2 3 4; do
    for thermal_conductivity in 1e-2 1e-3; do 
      echo "subsection Discretization" > current.prm
      echo "  set Temperature polynomial degree = $temperature_degree" >> current.prm
      echo "end" >> current.prm

      echo "subsection Material model" >> current.prm
      echo "  subsection Simpler model" >> current.prm
      echo "    set Thermal conductivity = $thermal_conductivity" >> current.prm
      echo "  end" >> current.prm
      echo "end" >> current.prm

      echo "subsection Mesh refinement" >> current.prm
      echo "  set Initial global refinement = $refinement" >> current.prm
      echo "end" >> current.prm


      echo "set Output directory = Q${temperature_degree}_refinement${refinement}_${thermal_conductivity}" >> current.prm
      echo "Starting Q${temperature_degree}_refinement${refinement}_${thermal_conductivity}"
      cat advection_in_annulus.prm current.prm | ./aspect --
    done
  done
done
