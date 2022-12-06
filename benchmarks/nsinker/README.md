The NSinker benchmark
=====================

This is the "NSinker" benchmark defined in

> Johann Rudi, Georg Stadler, Omar Ghattas.
> Weighted BFBT Preconditioner for Stokes Flow Problems with Highly
> Heterogeneous Viscosity
> SIAM J. Sci. Comput., 39(5), S272–S297.
> https://doi.org/10.1137/16M108450X

It creates a number of spherical high-viscosity, high-density sinking
spheres in a box geometry that provide a challenge for the Stokes
preconditioner. The difficulty of the problem is determined by the
number of sinkers and the viscosity contrast between sinkers and
background, which can be adjusted in the parameter file.
