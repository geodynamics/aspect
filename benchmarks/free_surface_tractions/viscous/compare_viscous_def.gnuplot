# Plot the maximum surface depression through time, due to the instantaneous
# cylindrical surface loading of a viscous medium. Compares ASPECT outputs
# with analytical solution.

set title "Depression of surface in response to loading (instaneous at t=0), over viscous medium"
set xlabel "Time [years]"
set ylabel "Maximum depression of surface [m]"

plot 'output_free_surface_viscous_cylinder_2D_loading/statistics' using 2:27 title 'Numerical (2-D)' with linesp
replot 'soln_z0.txt' using 1:2 title 'Analytical (Haskell, 1935a)'
