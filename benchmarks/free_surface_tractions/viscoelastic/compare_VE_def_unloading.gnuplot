# Plot the maximum surface depression through time, due to the linear
# cylindrical surface unloading of a viscoelastic medium. Compares ASPECT outputs
# with analytical solution.

set title "Depression of surface in response to linear unloading (from 0 to 1000 years), over viscoelastic medium"
set xlabel "Time [years]"
set ylabel "Maximum depression of surface [m]"

plot 'output_free_surface_VE_cylinder_2D_loading_unloading/statistics' using 2:28 title 'Numerical, no stress averaging' with linesp
replot 'soln2_z0.txt' using 1:2 title 'Analytical (Nakiboglu and Lambeck, 1982)'
