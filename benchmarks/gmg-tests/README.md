```{tags}
category:benchmark
feature:solver-comparison
```

(sec:benchmarks:gmg-tests)=
# Linear solver convergence tests

A collection of models to examine the performance and capabilities of
the GMG linear solver. Currently, the folder only contains a single test
set. This will be expanded in the future to include other geometries
and problems.


## Test 1: nsinker-box-3d

The 3d box (nsinker) problem with various combinations of boundary conditions:

- no-slip: all 6 sides no slip
- free-slip: all sides free slip
- partial-free-slip: left and front free slip, rest no slip
- partial-set: left: free slip; front: x set to smooth function, y=0; other: no slip
- open: no slip on all sides except top, which is open
- periodic: periodic in x direction; no slip otherwise
- free-surface: free surface on top with smooth topography, no slip bottom, free slip rest

The common test problem settings are:
- Q2-Q1
- tolerance: 1e-8
- 4 sinkers
- viscosity ratio: 1e4

Stokes GMRES iterations on the final mesh with 7.8m Stokes DoFs (10.3m total DoFs) with adaptive refinement:

| boundary condition | AMG | LS GMG | GC GMG | GC GMG (+CG AMG) |
|---|---|---|---|---|
| no-slip | 70 | 39 | 37 | 36 |
| free-slip | 172 | 73 | 61 | 36 |
| partial-free-slip | 105 | 53 | 51 | 45 |
| partial-set | 62 | 33 | 32 | 28 |
| open | 88 | 44 | 42 | 41 |
| periodic | 85 | - | 33 | 32 |
| free-surface | 195 | 78 | 67 | 39 |

GC GMG uses a Chebyshev coarse solver on the full hierarchy; GC GMG (+CG AMG)
uses AMG-preconditioned CG on a coarse mesh with at least 256 cells.

Currently unsupported cases:
- GMG-LS: no support for periodic boundaries
