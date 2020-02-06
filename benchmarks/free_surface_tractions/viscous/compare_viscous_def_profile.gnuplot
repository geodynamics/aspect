# Plot a radial profile of the surface depression at different timesteps, due 
# to the instantaneous cylindrical surface loading of a viscous medium.
# Compares ASPECT outputs with analytical solution.

set title "Depression of surface in response to loading (instaneous at t=0 years), over viscous medium"
set xlabel "Horizontal distance from load center, r [m]"
set ylabel "Depression of surface [m]"

plot for [i=2:8:2] 'output_free_surface_viscous_cylinder_2D_loading/topography.000'.i.'0' using 1:3 title 'Numerical (2-D), t='.i.'00 yrs'  
replot for [i=10:14:2] 'output_free_surface_viscous_cylinder_2D_loading/topography.00'.i.'0' using 1:3 title 'Numerical (2-D), t='.i.'00 yrs'

replot for [i=2:8:2] 'output_free_surface_viscous_cylinder_3D_loading/topography.000'.i.'0' using 1:3 title 'Numerical (3-D), t='.i.'00 yrs'
replot for [i=10:14:2] 'output_free_surface_viscous_cylinder_3D_loading/topography.00'.i.'0' using 1:3 title 'Numerical (3-D), t='.i.'00 yrs'

replot for [i=2:15:2] 'soln_'.i.'.txt' using 1:2 title 'Analytical, t='.i.'00 yrs' with lines
