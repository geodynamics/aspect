The NSinker benchmark
=====================

This is the "NSinker" benchmark defined in

> Weighted BFBT Preconditioner for Stokes Flow Problems with Highly Heterogeneous Viscosity,
> Johann Rudi, Georg Stadler, Omar Ghattas,
> SIAM Journal of Scientific Computing, Vol 39, Issue 5, 2017.
> https://doi.org/10.1137/16M108450X

It creates a number of spherical high-viscosity, high-density sinking
spheres in a box geometry that provide a challenge for the Stokes
preconditioner. The difficulty of the problem is determined by the
number of sinkers and the viscosity contrast between sinkers and
background, which can be adjusted in the parameter file.

Files for this benchmark are located in
[this directory](https://github.com/geodynamics/aspect/tree/main/benchmarks/nsinker).
