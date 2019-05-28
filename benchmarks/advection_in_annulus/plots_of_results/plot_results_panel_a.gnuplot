# Plot results from Advection  benchmark.
# Includes curves for average and maximum temperature, as well
# as an approximate analytical solution.
# Use ':%s/arg1/arg2/g' to replace all instances of arg1
# with arg2.

plot "results_Q2_beta0.0_cR0.0.dat" using 1:2  with linespoints, \
     "results_Q2_beta0.078_cR0.33.dat" using 1:2 with linespoints, \
     "results_Q2_beta0.052_cR0.11.dat" using 1:2  with linespoints, \
     "results_Q2_beta0.026_cR0.06.dat" using 1:2 with linespoints

set title "Convergence of Stabilization Parameters" font ',16'
set xlabel "Mesh Refinement" font ',12'
set ylabel "Average Temperature" font ',12'

pause -1

