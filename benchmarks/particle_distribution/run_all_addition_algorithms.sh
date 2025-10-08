#!/bin/bash

# Run all models in this directory
# > >> are used to append, cat concatenates


processes=4
ASPECT_EXEC="../../build/aspect"


# ------------------------------------Oscillating velocity------------------------------------ #

# Random addition, no point density function used
echo "subsection Prescribed Stokes solution" > current.prm
echo "set Model name = function" >> current.prm
echo "subsection Velocity function" >> current.prm
echo "  set Variable names      = x,y,t" >> current.prm
echo "  set Function constants  = velConstant=-0" >> current.prm
echo "  set Function expression = 0; (-0.5*sin(pi*t)) +velConstant" >> current.prm
echo " end" >> current.prm
echo "end" >> current.prm
echo "subsection Particles" >> current.prm
echo "  set Particle addition algorithm = random" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_addition/output-random" >> current.prm
cat addition_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --

# Random addition, no point density function used
echo "subsection Prescribed Stokes solution" > current.prm
echo "set Model name = function" >> current.prm
echo "subsection Velocity function" >> current.prm
echo "  set Variable names      = x,y,t" >> current.prm
echo "  set Function constants  = velConstant=-0" >> current.prm
echo "  set Function expression = 0; (-0.5*sin(pi*t)) +velConstant" >> current.prm
echo " end" >> current.prm
echo "end" >> current.prm
echo "subsection Particles" >> current.prm
echo "  set Particle addition algorithm = point density function" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_addition/output-point-density-function" >> current.prm
cat addition_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --
