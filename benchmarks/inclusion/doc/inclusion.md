(sec:benchmarks:inclusion)=
# The &ldquo;inclusion&rdquo; Stokes benchmark

The &ldquo;inclusion&rdquo; benchmark again solves a problem with a
discontinuous viscosity, but this time the viscosity is chosen in such a way
that the discontinuity is along a circle. This ensures that, unlike in the
[SolCx benchmark](sec:benchmarks:solcx), the discontinuity in the viscosity never
aligns to cell boundaries, leading to much larger difficulties in obtaining an
accurate representation of the pressure. Specifically, the almost
discontinuous pressure along this interface leads to oscillations in the
numerical solution. This can be seen in the visualizations shown in
{numref}`fig:inclusion1` and {numref}`fig:inclusion2`. As before, for details we refer to {cite:t}`duretz:etal:2011`. The
analytic solution against which we compare is given in {cite:t}`schmid:podladchikov:2003`. An extensive discussion of convergence properties is given in {cite:t}`kronbichler:etal:2012`.

```{figure-md} fig:inclusion1
<img src="inclusion-solution.*" />

The viscosity field when interpolated onto the mesh (internally, the exact viscosity field is used): large inside a circle and small outside it. Overlaid to it a vector field showing the velocity.
```

```{figure-md} fig:inclusion2
<img src="inclusion-solution-pressure.*" />

The pressure with its oscillations along the interface. The oscillations become more localized as the mesh is refined.
```

The benchmark can be run using the parameter files in
[benchmarks/inclusion/](https://github.com/geodynamics/aspect/blob/main/benchmarks/inclusion). The material model, boundary condition, and
postprocessor are defined in [benchmarks/inclusion/inclusion.cc](https://github.com/geodynamics/aspect/blob/main/benchmarks/inclusion/inclusion.cc).
Consequently, this code needs to be compiled into a shared lib before you can
run the tests.


``` {literalinclude} ./inclusion.prm
```
