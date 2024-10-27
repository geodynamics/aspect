The Spherical Shell NSinker benchmark
=====================================

[This folder](https://github.com/geodynamics/aspect/tree/main/benchmarks/nsinker_spherical_shell)
contains the "Spherical Shell NSinker" benchmark, which is
loosely based on the NSinker benchmark defined in

> Weighted BFBT Preconditioner for Stokes Flow Problems with Highly Heterogeneous Viscosity,
> Johann Rudi, Georg Stadler, Omar Ghattas,
> SIAM Journal of Scientific Computing, Vol 39, Issue 5, 2017.
> https://doi.org/10.1137/16M108450X

but adapted to a spherical shell geometry. See the "nsinker" benchmark
in this repository for the version in a unit box. The motivation for
creating this benchmark was the desire to test linear solver
performance with a spherical geometry, which is often used in practice
and that turns out to be more challenging.

This benchmark is defined and used in

> Efficient distributed matrix-free multigrid methods on locally refined meshes for FEM computations,
> Peter Munch, Timo Heister, Laura Prieto Saavedra, Martin Kronbichler,
> ACM Transactions on Parallel Computing, 2023, 10, 1-38.
> https://doi.org/10.1145/3580314

We consider a spherical shell with inner radius 0.54 and outer radius
1.0 with N spherical, high-viscosity, high-density sinking spheres at
random (but predetermined in the code) locations in the
domain. Additionally, a transition zone with distance 0.89583 from the
origin (or about 660km depth if this is assumed to be a
non-dimensionalized Earth-like computation) with a sharp viscosity
jump can be enabled.
