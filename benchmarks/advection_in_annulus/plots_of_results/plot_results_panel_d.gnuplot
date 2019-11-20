# Plot results from Advection  benchmark.
# Includes curves for average and maximum temperature, as well
# as an approximate analytical solution.
# Use ':%s/arg1/arg2/g' to replace all instances of arg1
# with arg2.

plot "results_Q1_refinement4_beta0.052.dat" using 1:3  with linespoints, \
     "results_Q2_refinement4_beta0.052.dat" using 1:3 with linespoints, \
     "results_Q3_refinement4_beta0.052.dat" using 1:3 with linespoints

set title "Convergence of Temperature Polynomials due to cR" font ',16'
set xlabel "cR" font ',12'
set ylabel "Maximum Temperature" font ',12'

pause -1

