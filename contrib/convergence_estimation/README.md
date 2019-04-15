# Convergence Estimation

This is a python script intended to work with the `state quadrature`
post-processing plugin to allow comparison of solutions using `$L_p$` norms for
purposes such as estimating convergence rate of problems that do not have an
exact solution.

## Usage

### Caveats

Due to MPI mesh divisions and other factors, even if the originating meshes are
identical, the order of the points cannot be assumed to be consistent between
files. Therefore, the order must be modified by the postprocessor to achieve
a consistent order. Doing this by a comparison method would be `$O(n^2)$` in
the number of points assuming that a matching point can be identified using
only the two points in question. Given that the post-processor is likely to be
applied to large data-sets a faster method would be preferable.

The current script has implemented a quicker algorithm to obtain a consistent
ordering for geometries based on the Box `GeometryModel`, but is not functional
for general geometries as it relies on the structure of the mesh to produce
consistent results.
