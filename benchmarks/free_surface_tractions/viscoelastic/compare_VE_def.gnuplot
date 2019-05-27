# Plot the maximum surface depression through time, due to the instantaneous
# cylindrical surface loading of a viscoelastic medium. Compares ASPECT outputs
# with analytical solution.

set title "Depression of surface in response to loading (inst. from 0 to 1000 years), over viscoelastic medium"
set xlabel "Time [years]"
set ylabel "Maximum depression of surface [m]"

plot 'output_free_surface_VE_cylinder_2D_loading/statistics' using 2:27 title 'Numerical, no stress averaging' with linesp
replot 'output_free_surface_VE_cylinder_2D_loading_fixed_elastic_dt/statistics' using 2:27 title 'Numerical, with stress averaging' with linesp
replot 'soln_z0.txt' using 1:2 title 'Analytical (Nakiboglu and Lambeck, 1982)'
