```{tags}
category:benchmark
feature:solver-comparison
```

(sec:benchmarks:gmg-tests)=
# GMG linear solver convergence tests

## Test 1: box-3d-bc/

3d box (nsinker) with various combinations of boundary conditions:

no-slip: all 6 sides no slip
free-slip: all sides free slip
partial-free-slip: left and front free slip, rest no slip
partial-set: left: free slip; front: x set to smooth function, y=0; other: no slip
open: no slip on all sides except top, which is open
periodic: periodic in x direction; no slip otherwise
free-surface: free surface on top with smooth topography, no slip bottom, free slip rest

Stokes: Q2Q1, 7.8m Stokes DoFs (10.3m total DoFs) on the final mesh (adaptively refined)
tolerance: 1e-8
4 sinkers, viscosity ratio: 1e4

Stokes GMRES iterations on the final (largest) mesh:

| boundary condition | GMG global coarsening | GMG local smoothing | AMG |
|---|---|---|---|
| no-slip | 37 | 39 | 70 |
| free-slip | 61 | 73 | 172 |
| partial-free-slip | 51 | 53 | 105 |
| partial-set | 32 | 33 | 62 |
| open | 42 | 44 | 88 |
| periodic | 33 | - | 85 |
| free-surface | - | 78 | 195 |

Note:
- GMG-LS: no support for periodic boundaries
- GMG-GC: No support for mesh deformation (free surface)
