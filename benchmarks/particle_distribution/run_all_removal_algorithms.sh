#!/bin/bash

# Run all removal benchmarks in this directory
# > >> are used to append, cat concatenates


processes=4
ASPECT_EXEC="../../build/aspect"

# ------------------------------------Oscillating velocity------------------------------------ #

# Random removal, no point density function used
echo "subsection Prescribed Stokes solution" > current.prm
echo "set Model name = function" >> current.prm
echo "subsection Velocity function" >> current.prm
echo "  set Variable names      = x,y,t" >> current.prm
echo "  set Function expression = 0; (-0.5*sin(pi*t))" >> current.prm
echo " end" >> current.prm
echo "end" >> current.prm
echo "subsection Particles" >> current.prm
echo "  set Particle removal algorithm = random" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_removal/output-random" >> current.prm
cat removal_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --

# Cutoff_c1 (the default kernel function)
echo "subsection Prescribed Stokes solution" > current.prm
echo "set Model name = function" >> current.prm
echo "subsection Velocity function" >> current.prm
echo "  set Variable names      = x,y,t" >> current.prm
echo "  set Function expression = 0; (-0.5*sin(pi*t))" >> current.prm
echo " end" >> current.prm
echo "end" >> current.prm
echo "subsection Particles" >> current.prm
echo "  set Particle removal algorithm = point density function" >> current.prm
echo "  set Point density kernel function = cutoff c1 dealii" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_removal/output-cutoff-c1" >> current.prm
cat removal_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --

# Cutoff_w1
echo "subsection Prescribed Stokes solution" > current.prm
echo "set Model name = function" >> current.prm
echo "subsection Velocity function" >> current.prm
echo "  set Variable names      = x,y,t" >> current.prm
echo "  set Function expression = 0; (-0.5*sin(pi*t))" >> current.prm
echo " end" >> current.prm
echo "end" >> current.prm
echo "subsection Particles" >> current.prm
echo "  set Particle removal algorithm = point density function" >> current.prm
echo "  set Point density kernel function = cutoff w1 dealii" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_removal/output-cutoff-w1" >> current.prm
cat removal_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --


# Uniform kernel function
echo "subsection Prescribed Stokes solution" > current.prm
echo "set Model name = function" >> current.prm
echo "subsection Velocity function" >> current.prm
echo "  set Variable names      = x,y,t" >> current.prm
echo "  set Function expression = 0; (-0.5*sin(pi*t))" >> current.prm
echo " end" >> current.prm
echo "end" >> current.prm
echo "subsection Particles" >> current.prm
echo "  set Particle removal algorithm = point density function" >> current.prm
echo "  set Point density kernel function = uniform" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_removal/output-uniform" >> current.prm
cat removal_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --


# Gaussian kernel function
echo "subsection Prescribed Stokes solution" > current.prm
echo "set Model name = function" >> current.prm
echo "subsection Velocity function" >> current.prm
echo "  set Variable names      = x,y,t" >> current.prm
echo "  set Function expression = 0; (-0.5*sin(pi*t))" >> current.prm
echo " end" >> current.prm
echo "end" >> current.prm
echo "subsection Particles" >> current.prm
echo "  set Particle removal algorithm = point density function" >> current.prm
echo "  set Point density kernel function = gaussian" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_removal/output-gaussian" >> current.prm
cat removal_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --


# Triangular kernel function
echo "subsection Prescribed Stokes solution" > current.prm
echo "set Model name = function" >> current.prm
echo "subsection Velocity function" >> current.prm
echo "  set Variable names      = x,y,t" >> current.prm
echo "  set Function expression = 0; (-0.5*sin(pi*t))" >> current.prm
echo " end" >> current.prm
echo "end" >> current.prm
echo "subsection Particles" >> current.prm
echo "  set Particle removal algorithm = point density function" >> current.prm
echo "  set Point density kernel function = triangular" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_removal/output-triangular" >> current.prm
cat removal_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --

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
echo "  set Particle removal algorithm = random" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_removal/output-random-constant-velocity" >> current.prm
cat removal_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --

# Cutoff-c1 kernel function
echo "subsection Prescribed Stokes solution" > current.prm
echo "set Model name = function" >> current.prm
echo "subsection Velocity function" >> current.prm
echo "  set Variable names      = x,y" >> current.prm
echo "  set Function constants  = velSlow=-0.1" >> current.prm
echo "  set Function expression = 0; velSlow" >> current.prm
echo " end" >> current.prm
echo "end" >> current.prm
echo "subsection Particles" >> current.prm
echo "  set Particle removal algorithm = point density function" >> current.prm
echo "  set Point density kernel function = cutoff c1 dealii" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_removal/output-cutoff-c1-constant-velocity" >> current.prm
cat removal_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --

# Cutoff-w1 kernel function
echo "subsection Prescribed Stokes solution" > current.prm
echo "set Model name = function" >> current.prm
echo "subsection Velocity function" >> current.prm
echo "  set Variable names      = x,y" >> current.prm
echo "  set Function constants  = velSlow=-0.1" >> current.prm
echo "  set Function expression = 0; velSlow" >> current.prm
echo " end" >> current.prm
echo "end" >> current.prm
echo "subsection Particles" >> current.prm
echo "  set Particle removal algorithm = point density function" >> current.prm
echo "  set Point density kernel function = cutoff w1 dealii" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_removal/output-cutoff-w1-constant-velocity" >> current.prm
cat removal_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --


# Uniform kernel function
echo "subsection Prescribed Stokes solution" > current.prm
echo "set Model name = function" >> current.prm
echo "subsection Velocity function" >> current.prm
echo "  set Variable names      = x,y" >> current.prm
echo "  set Function constants  = velSlow=-0.1" >> current.prm
echo "  set Function expression = 0; velSlow" >> current.prm
echo " end" >> current.prm
echo "end" >> current.prm
echo "subsection Particles" >> current.prm
echo "  set Particle removal algorithm = point density function" >> current.prm
echo "  set Point density kernel function = uniform" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_removal/output-uniform-constant-velocity" >> current.prm
cat removal_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --


# Gaussian kernel function
echo "subsection Prescribed Stokes solution" > current.prm
echo "set Model name = function" >> current.prm
echo "subsection Velocity function" >> current.prm
echo "  set Variable names      = x,y" >> current.prm
echo "  set Function constants  = velSlow=-0.1" >> current.prm
echo "  set Function expression = 0; velSlow" >> current.prm
echo " end" >> current.prm
echo "end" >> current.prm
echo "subsection Particles" >> current.prm
echo "  set Particle removal algorithm = point density function" >> current.prm
echo "  set Point density kernel function = gaussian" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_removal/output-gaussian-constant-velocity" >> current.prm
cat removal_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --


# Triangular kernel function
echo "subsection Prescribed Stokes solution" > current.prm
echo "set Model name = function" >> current.prm
echo "subsection Velocity function" >> current.prm
echo "  set Variable names      = x,y" >> current.prm
echo "  set Function constants  = velSlow=-0.1" >> current.prm
echo "  set Function expression = 0; velSlow" >> current.prm
echo " end" >> current.prm
echo "end" >> current.prm
echo "subsection Particles" >> current.prm
echo "  set Particle removal algorithm = point density function" >> current.prm
echo "  set Point density kernel function = triangular" >> current.prm
echo "end" >> current.prm
echo "set Output directory = output_removal/output-triangular-constant-velocity" >> current.prm
cat removal_algorithm_benchmarks.prm current.prm | mpirun -np $processes $ASPECT_EXEC --

# Remove the temporary .prm file
rm current.prm

# ------------------------------------Calculate Averages--------------------------------------- #

python3 calculate_removal_algorithms_statistics.py

# ----------------------------------------Visual Plots----------------------------------------- #

gnuplot plot_all_removal_algorithms.gnuplot
