```{tags}
category:benchmark
```

(sec:benchmarks:viscous_half_space_loading)=
# Viscous half-space loading
Benchmark of infinite viscous half-space loaded/unloaded by an
 axisymmetric cylinder, using the solution of: Haskell, N.A.,
 The Motion of a Viscous Fluid Under a Surface Load (1935), Physics
 6, 265; doi:10.1063/1.1745329

This benchmark is identical to the viscoelastic benchmarks, in the
 limit of negligible elasticity.

The benchmark compares his solution for the surface displacement
 following the application of an instantaneous cylindrical surface
 load (Eq. 3.34) with the deformation of the free surface in ASPECT
 after imposing surface boundary tractions.

A surface pressure of rho_l*g*H0 (where rho_l is the load density,
 H0 is the load height) is applied instantaneously on the surface
 for r<r0 (where r0 is the load radius), for t>0. This is done both
 in a 2-D and 3-D geometry (by symmetry the load is centered on the
 left boundary or left/front corner). The input files are:
   'free_surface_viscous_cylinder_2D_loading.prm'
   'free_surface_viscous_cylinder_3D_loading.prm'

Over long timescales, the surface deformation beneath the load (i.e.
 ice sheet) should converge to H0*rho_l/rho_s, where rho_s is the
 density of the medium (i.e. mantle). The deformation outside the
 load should converge to zero. The numerical runs have slightly
 smaller deformations due to the pressure from the far vertical (right)
 boundary, this effect is lessened by increasing the horizontal domain
 size.

The 'topography' output files may be compared against an analytical
 solution by running the script '../viscoelastic/compare_VE_soln.py',
 which outputs the viscous limit of the Nakiboglu and Lambeck (1982)
 solution for an instantaneous cylindrical load on an infinite half-
 space. These outputs ('soln_XX.txt') may be compared against the
 ASPECT 'topography' output files, using the provided gnuplot script
 ('compare_viscous_def.gnuplot' for maximum surface deflection
 through time, 'compare_viscous_def_profile.gnuplot' for deflection
 of profile through time).

Note that while the analytical and numerical results for the deflection
 of the surface agree well near the center of the load (left boundary),
 the solutions do not match as well on the right (free-slip) boundary.
 The numerical free surface rides up vertically against the right
 boundary to conserve mass, while the analytical solution assumes an
 infinite half-space, predicting near-zero displacement far from the
 load center. This is also true for the visoelastic benchmark. This problem
 is less pronounced in 3-D as the extra mass may be distributed over
 a larger area. An open far (right/back) boundary resolves this problem
 in 3-D.

The solutions match well in 3-D. In 2-D, the geometries of the loading
 function are different (Cartesian in ASPECT vs cylindrical analytically).
 As such, the agreement in 2-D breaks down for small r0 (load width) or
 if the right boundary is placed further away.
