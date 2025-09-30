#!/bin/bash

# Run all addition benchmarks in this directory
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

# Histogram addition with default granularity
echo "subsection Prescribed Stokes solution" > current.prm
echo "set Model name = function" >> current.prm
echo "subsection Velocity function" >> current.prm
echo "  set Variable names      = x,y,t" >> current.prm
echo "  set Function constants  = velConstant=-0" >> current.prm
echo "  set Function expression = 0; (-0.5*sin(pi*t)) +velConstant" >> current.prm
echo " end" >> current.prm
echo "end" >> current.prm
echo "subsection Particles" >> current.prm
echo "  set Particle addition algorithm = histogram" >> current.prm
echo "  set Addition granularity = 8" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_addition/output-histogram" >> current.prm
cat addition_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --

# Point density function, Cutoff_c1 (the default kernel function)
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
echo "set Output directory = output_addition/output-cutoff-c1" >> current.prm
cat addition_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --

# Point density function, Cutoff_w1
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
echo "  set Point density kernel function = cutoff w1 dealii" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_addition/output-cutoff-w1" >> current.prm
cat addition_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --

# Point density function, uniform
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
echo "  set Point density kernel function = uniform" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_addition/output-uniform" >> current.prm
cat addition_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --

# Point density function, gaussian
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
echo "  set Point density kernel function = gaussian" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_addition/output-gaussian" >> current.prm
cat addition_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --

# Point density function, triangular
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
echo "  set Point density kernel function = triangular" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_addition/output-triangular" >> current.prm
cat addition_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --

# ------------------------------------Constant velocity------------------------------------ #

# Random removal, no point density function used
echo "subsection Prescribed Stokes solution" > current.prm
echo "set Model name = function" >> current.prm
echo "subsection Velocity function" >> current.prm
echo "  set Variable names      = x,y" >> current.prm
echo "  set Function constants  = velSlow=-0.1" >> current.prm
echo "  set Function expression = 0; velSlow" >> current.prm
echo " end" >> current.prm
echo "end" >> current.prm
echo "subsection Particles" >> current.prm
echo "  set Particle addition algorithm = random" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_addition/output-random-constant-velocity" >> current.prm
cat addition_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --

# Histogram addition with default granularity
echo "subsection Prescribed Stokes solution" > current.prm
echo "set Model name = function" >> current.prm
echo "subsection Velocity function" >> current.prm
echo "  set Variable names      = x,y" >> current.prm
echo "  set Function constants  = velSlow=-0.1" >> current.prm
echo "  set Function expression = 0; velSlow" >> current.prm
echo " end" >> current.prm
echo "end" >> current.prm
echo "subsection Particles" >> current.prm
echo "  set Particle addition algorithm = histogram" >> current.prm
echo "  set Addition granularity = 8" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_addition/output-histogram-constant-velocity" >> current.prm
cat addition_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --

# Point density function, Cutoff_c1 (the default kernel function)
echo "subsection Prescribed Stokes solution" > current.prm
echo "set Model name = function" >> current.prm
echo "subsection Velocity function" >> current.prm
echo "  set Variable names      = x,y" >> current.prm
echo "  set Function constants  = velSlow=-0.1" >> current.prm
echo "  set Function expression = 0; velSlow" >> current.prm
echo " end" >> current.prm
echo "end" >> current.prm
echo "subsection Particles" >> current.prm
echo "  set Particle addition algorithm = point density function" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_addition/output-cutoff-c1-constant-velocity" >> current.prm
cat addition_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --

# Point density function, Cutoff_w1 (the default kernel function)
echo "subsection Prescribed Stokes solution" > current.prm
echo "set Model name = function" >> current.prm
echo "subsection Velocity function" >> current.prm
echo "  set Variable names      = x,y" >> current.prm
echo "  set Function constants  = velSlow=-0.1" >> current.prm
echo "  set Function expression = 0; velSlow" >> current.prm
echo " end" >> current.prm
echo "end" >> current.prm
echo "subsection Particles" >> current.prm
echo "  set Particle addition algorithm = point density function" >> current.prm
echo "  set Point density kernel function = cutoff w1 dealii" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_addition/output-cutoff-w1-constant-velocity" >> current.prm
cat addition_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --

# Point density function, uniform
echo "subsection Prescribed Stokes solution" > current.prm
echo "set Model name = function" >> current.prm
echo "subsection Velocity function" >> current.prm
echo "  set Variable names      = x,y" >> current.prm
echo "  set Function constants  = velSlow=-0.1" >> current.prm
echo "  set Function expression = 0; velSlow" >> current.prm
echo " end" >> current.prm
echo "end" >> current.prm
echo "subsection Particles" >> current.prm
echo "  set Particle addition algorithm = point density function" >> current.prm
echo "  set Point density kernel function = uniform" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_addition/output-uniform-constant-velocity" >> current.prm
cat addition_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --

# Point density function, gaussian
echo "subsection Prescribed Stokes solution" > current.prm
echo "set Model name = function" >> current.prm
echo "subsection Velocity function" >> current.prm
echo "  set Variable names      = x,y" >> current.prm
echo "  set Function constants  = velSlow=-0.1" >> current.prm
echo "  set Function expression = 0; velSlow" >> current.prm
echo " end" >> current.prm
echo "end" >> current.prm
echo "subsection Particles" >> current.prm
echo "  set Particle addition algorithm = point density function" >> current.prm
echo "  set Point density kernel function = gaussian" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_addition/output-gaussian-constant-velocity" >> current.prm
cat addition_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --

# Point density function, triangular
echo "subsection Prescribed Stokes solution" > current.prm
echo "set Model name = function" >> current.prm
echo "subsection Velocity function" >> current.prm
echo "  set Variable names      = x,y" >> current.prm
echo "  set Function constants  = velSlow=-0.1" >> current.prm
echo "  set Function expression = 0; velSlow" >> current.prm
echo " end" >> current.prm
echo "end" >> current.prm
echo "subsection Particles" >> current.prm
echo "  set Particle addition algorithm = point density function" >> current.prm
echo "  set Point density kernel function = triangular" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_addition/output-triangular-constant-velocity" >> current.prm
cat addition_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --

# Remove the temporary .prm file
rm current.prm

# ------------------------------------Calculate Averages--------------------------------------- #

python3 calculate_addition_algorithms_statistics.py

# ----------------------------------------Visual Plots----------------------------------------- #

gnuplot plot_all_addition_algorithms.gnuplot
