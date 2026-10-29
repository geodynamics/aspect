# Dynamic Topography Adjoint Benchmark: 2D Spherical Shell

This benchmark computes cellwise density and viscosity kernels for the
dynamic-topography objective in an instantaneous spherical-shell model.
Selected physical-property cell perturbations are validated with centered
finite differences.

Run the benchmark from a build directory with:

```bash
./aspect ../benchmarks/adjoint_dynamic_topography_2d_spherical/adjoint_dynamic_topography.prm
```

Then create the cellwise density/viscosity kernel plot:

```bash
../benchmarks/adjoint_dynamic_topography_2d_spherical/scripts/plot_cellwise_kernels.py \
  output-adjoint-dynamic-topography-2d-spherical \
  --output ../benchmarks/adjoint_dynamic_topography_2d_spherical/doc/cellwise_density_viscosity_kernels.png
```

To run the finite-difference validation, first build the benchmark material
wrapper:

```bash
cmake -S ../benchmarks/adjoint_dynamic_topography_2d_spherical \
  -B adjoint_dynamic_topography_2d_spherical_plugin \
  -D Aspect_DIR=$(pwd)
cmake --build adjoint_dynamic_topography_2d_spherical_plugin
cp adjoint_dynamic_topography_2d_spherical_plugin/libadjoint_dynamic_topography_2d_spherical.release.so \
  ../benchmarks/adjoint_dynamic_topography_2d_spherical/
```

Then run:

```bash
../benchmarks/adjoint_dynamic_topography_2d_spherical/scripts/run_fd_validation.py \
  --aspect ./aspect \
  --selection largest \
  --cells-per-property 32
```

The validation wrapper perturbs physical density or viscosity directly in a
selected active cell after material evaluation. Therefore the finite-difference
directional derivative is compared with the corresponding physical-property
adjoint kernel times cell volume.
