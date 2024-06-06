# The Rayleigh-Taylor instability with a free surface

*This section was contributed by Anne Glerum, Ian Rose and Timo Heister.*

The prm file kaus_rayleigh_taylor_instability.prm runs one of four cases
to test the free surface stabilization as described in Kaus et al. (2010).
1. No stabilization, fixed timestep of 2500 yr, vertical projection.
2. No stabilization, fixed timestep of 10000 yr, vertical projection.
3. With stabilization, fixed timestep of 10000 yr, vertical projection.
4. With stabilization, fixed timestep of 10000 yr, normal projection.

Plotting the maximum topography with plot_topography.gnuplot (e.g., {numref}`fig:topography`)
and the maximum depth of the drip (e.g., {numref}`fig:topography`) along the left vertical
boundary with plot_drip_depth.gnuplot shows that:
- Case 2 leads to 'drunken sailor' behavior.
- Case 3 is stable again and matches the results of Kaus et al. (2010) well.
- The topography between case 3 and 4 differs substantially.

The prm file rose_rayleigh_taylor_instability.prm is based on Rose et al. (2017)
and shows how adaptive mesh refinement can be used to resolve the layer interface
and the surface at approximately 1 km resolution instead of the global resolution
of ~8 km in the first prm. Plotting the topography with plot_topography.gnuplot
shows the effect of this higher resolution.

```{figure-md} fig:drip_depth
<img src="drip_depth.*" />

Drip depth of an ASPECT simulation with free surface stabilization, normal projection
of the free surface and a maximum timestep of 5000 yr. The drip depth matches well with
the results published in Kaus et al. (2010).
```

```{figure-md} fig:topography
<img src="topography.*" />

Maximum topography of an ASPECT simulation with free surface stabilization, normal projection
of the free surface and a maximum timestep of 5000 yr.
```
