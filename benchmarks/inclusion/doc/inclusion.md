#### The "inclusion" Stokes benchmark

The "inclusion" benchmark again solves a problem with a
discontinuous viscosity, but this time the viscosity is chosen in such a way
that the discontinuity is along a circle. This ensures that, unlike in the
SolCx benchmark discussed above, the discontinuity in the viscosity never
aligns to cell boundaries, leading to much larger difficulties in obtaining an
accurate representation of the pressure. Specifically, the almost
discontinuous pressure along this interface leads to oscillations in the
numerical solution. This can be seen in the visualizations shown in
Fig. [2]. As before, for details we refer to (Duretz et al. 2011). The
analytic solution against which we compare is given in (Schmid and
Podladchikov 2003). An extensive discussion of convergence properties is given
in (Kronbichler, Heister, and Bangerth 2012).

<div class="center">

<img src="cookbooks/benchmarks/inclusion/doc/inclusion-solution.png" title="fig:" id="fig:inclusion" style="width:45.0%" alt="Inclusion Stokes benchmark. Left: The viscosity field when interpolated onto the mesh (internally, the &#x201C;exact&#x201D; viscosity field &#x2013; large inside a circle, small outside &#x2013; is used), and overlaid to it some velocity vectors. Right: The pressure with its oscillations along the interface. The oscillations become more localized as the mesh is refined." />
<img src="cookbooks/benchmarks/inclusion/doc/inclusion-solution-pressure.png" title="fig:" id="fig:inclusion" style="width:45.0%" alt="Inclusion Stokes benchmark. Left: The viscosity field when interpolated onto the mesh (internally, the &#x201C;exact&#x201D; viscosity field &#x2013; large inside a circle, small outside &#x2013; is used), and overlaid to it some velocity vectors. Right: The pressure with its oscillations along the interface. The oscillations become more localized as the mesh is refined." />

</div>

The benchmark can be run using the parameter files in
[benchmarks/inclusion/]. The material model, boundary condition, and
postprocessor are defined in [benchmarks/inclusion/inclusion.cc].
Consequently, this code needs to be compiled into a shared lib before you can
run the tests.

``` prmfile
```

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-DMGT11" class="csl-entry">

Duretz, T., D. A. May, T. V. Garya, and P. J. Tackley. 2011.
"Discretization Errors and Free Surface Stabilization in the Finite
Difference and Marker-in-Cell Method for Applied Geodynamics: A Numerical
Study." *Geoch. Geoph. Geosystems* 12: Q07004/1&ndash;26.

</div>

<div id="ref-KHB12" class="csl-entry">

Kronbichler, M., T. Heister, and W. Bangerth. 2012. "High Accuracy
Mantle Convection Simulation Through Modern Numerical Methods."
*Geophysical Journal International* 191: 12&ndash;29.
<https://doi.org/10.1111/j.1365-246X.2012.05609.x>.

</div>

<div id="ref-SP03" class="csl-entry">

Schmid, D. W., and Y. Y. Podladchikov. 2003. "Analytical Solutions for
Deformable Elliptical Inclusions in General Shear." *Geophysical Journal
International* 155 (1): 269&ndash;88.

</div>

</div>

  [2]: #fig:inclusion
  [benchmarks/inclusion/]: benchmarks/inclusion/
  [benchmarks/inclusion/inclusion.cc]: benchmarks/inclusion/inclusion.cc
