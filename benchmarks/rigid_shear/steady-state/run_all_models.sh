#!/bin/bash

# Run all models in this directory. Note that the highest resolutions and ppd require
# a high-end laptop or better to run (32 GB RAM).

processes=8
ASPECT_EXEC="../plugin/aspect"

for stokes_degree in 2 3; do
  for discontinuous_pressure in false; do
    for refinement in 3 4 5 6; do
      for particles_per_direction in 4 5 6 7 10 15 20 32 45 64 80; do
        for generator in "reference cell"; do
        for interpolator in "cell average" "bilinear least squares"; do
        for integrator in "rk2" "rk4"; do

        echo "subsection Discretization" > current.prm
        echo "  set Stokes velocity polynomial degree = $stokes_degree" >> current.prm
        echo "  set Use locally conservative discretization = $discontinuous_pressure" >> current.prm
        echo "end" >> current.prm

        echo "subsection Particles" >> current.prm
        echo "  set Particle generator name = $generator" >> current.prm

        # If using longer runtimes it is necessary to add particles to
        # prevent empty cells. Compute an appropriate number.
        #ppc=`expr $particles_per_direction \* $particles_per_direction`
        #min_ppc=`expr $ppc / 2`
        #echo "  set Minimum particles per cell = $min_ppc" >> current.prm
        #echo "  set Load balancing strategy = add particles" >> current.prm

        echo "  set Interpolation scheme = $interpolator" >> current.prm
        echo "  set Integration scheme = $integrator" >> current.prm
        echo "  subsection Generator" >> current.prm
        echo "    subsection Reference cell" >> current.prm
        echo "      set Number of particles per cell per direction = $particles_per_direction" >> current.prm
        echo "    end" >> current.prm
        echo "  end" >> current.prm
        echo "end" >> current.prm

        echo "subsection Mesh refinement" >> current.prm
        echo "  set Initial global refinement = $refinement" >> current.prm
        echo "end" >> current.prm


        echo "set Output directory = Q${stokes_degree}_P${discontinuous_pressure}_refinement${refinement}_${particles_per_direction}_${generator// /_}_${interpolator// /_}_${integrator}" >> current.prm
        echo "Starting Q${stokes_degree}_P${discontinuous_pressure}_refinement${refinement}_${particles_per_direction}_${integrator}"
        cat rigid_shear.prm current.prm | mpirun -np $processes $ASPECT_EXEC --

        done
        done
        done
      done
    done
  done
done
