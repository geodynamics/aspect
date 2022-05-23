#### The &ldquo;inclusion&rdquo; Stokes benchmark

The &ldquo;inclusion&rdquo; benchmark again solves a problem with a
discontinuous viscosity, but this time the viscosity is chosen in such a way
that the discontinuity is along a circle. This ensures that, unlike in the
SolCx benchmark discussed above, the discontinuity in the viscosity never
aligns to cell boundaries, leading to much larger difficulties in obtaining an
accurate representation of the pressure. Specifically, the almost
discontinuous pressure along this interface leads to oscillations in the
numerical solution. This can be seen in the visualizations shown in
{numref}`fig:inclusion1` and {numref}`fig:inclusion2`. As before, for details we refer to {cite}`DMGT11`. The
analytic solution against which we compare is given in {cite}`SP03`. An extensive discussion of convergence properties is given in {cite}`kronbichler:etal:2012`.

```{figure-md} fig:inclusion1
<img src="inclusion-solution.*" alt="The viscosity field when interpolated onto the mesh (internally, the &#x201C;exact&#x201D; viscosity field &#x2013; large inside a circle, small outside &#x2013; is used), and overlaid to it some velocity vectors." />

The viscosity field when interpolated onto the mesh (internally, the &#x201C;exact&#x201D; viscosity field &#x2013; large inside a circle, small outside &#x2013; is used), and overlaid to it some velocity vectors.
```
```{figure-md} fig:inclusion2
<img src="inclusion-solution-pressure.*"  alt="The pressure with its oscillations along the interface. The oscillations become more localized as the mesh is refined." />

The pressure with its oscillations along the interface. The oscillations become more localized as the mesh is refined.
```
The benchmark can be run using the parameter files in
[benchmarks/inclusion/](https://github.com/geodynamics/aspect/blob/main/benchmarks/inclusion). The material model, boundary condition, and
postprocessor are defined in [benchmarks/inclusion/inclusion.cc](https://github.com/geodynamics/aspect/blob/main/benchmarks/inclusion/inclusion.cc).
Consequently, this code needs to be compiled into a shared lib before you can
run the tests.


``` {literalinclude} ./inclusion.prm
```
