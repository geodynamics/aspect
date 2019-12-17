The prm file runs one of four cases to test the free surface stabilization
as described in Kaus et al. (2010) and Rose et al. (2017). 
1. No stabilization, fixed timestep of 2500 yr, vertical projection.
2. No stabilization, fixed timestep of 5000 yr, vertical projection.
3. With stabilization, fixed timestep of 5000 yr, vertical projection.
4. With stabilization, fixed timestep of 5000 yr, normal projection.

Plotting the maximum topography with plot_topography.gnuplot and the maximum
depth of the drip along the left vertical boundary with plot_drip_depth.gnuplot
shows that:
- Case 2 leads to 'drunken sailor' behavior
- Case 3 is stable again
- The topography between case 3 and 4 differs substantially.
