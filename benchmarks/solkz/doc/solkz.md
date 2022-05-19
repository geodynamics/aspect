#### The SolKz Stokes benchmark

The SolKz benchmark is another variation on the same theme as the SolCx
benchmark above: it solves a Stokes problem with a spatially variable
viscosity, but this time the viscosity is not a discontinuous function.
Instead, it grows exponentially with the vertical coordinate so that its
overall variation is again $10^6$. The forcing is again chosen by imposing a
spatially variable density variation. For details, refer again to (Duretz et
al. 2011).

The following input file, only a small variation of the one in the previous
section, solves the benchmark (see [benchmarks/solkz/][]):

``` prmfile
```

The output when running on this parameter file looks similar to the one shown
for the SolCx case. The solution when computed with one more level of global
refinement is visualized in Fig.&nbsp;[2][]. The velocity solution computed
with three more levels of global refinement and plotted over the viscosity
field is shown in Fig.&nbsp;[3][].

<div class="center">

<img src="cookbooks/benchmarks/solkz/doc/solkz-solution.png" title="fig:" id="fig:solkz" style="width:45.0%" alt="SolKz Stokes benchmark. Left: The density perturbation field overlaid with velocity vectors. The viscosity grows exponentially in the vertical direction, leading to small velocities at the top despite the large density variations. Right: The pressure." />
<img src="cookbooks/benchmarks/solkz/doc/solkz-solution-pressure.png" title="fig:" id="fig:solkz" style="width:45.0%" alt="SolKz Stokes benchmark. Left: The density perturbation field overlaid with velocity vectors. The viscosity grows exponentially in the vertical direction, leading to small velocities at the top despite the large density variations. Right: The pressure." />

</div>

<div class="center">

<figure>
<img src="cookbooks/benchmarks/solkz/doc/solkz-solution-viscosity.png" id="fig:solkz2" style="width:70.0%" alt="SolKz Stokes benchmark. Another view of the velocity vectors, this time plotted over the viscosity field." /><figcaption aria-hidden="true"><em>SolKz Stokes benchmark. Another view of the velocity vectors, this time plotted over the viscosity field.</em></figcaption>
</figure>

</div>

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-DMGT11" class="csl-entry">

Duretz, T., D. A. May, T. V. Garya, and P. J. Tackley. 2011.
&ldquo;Discretization Errors and Free Surface Stabilization in the Finite
Difference and Marker-in-Cell Method for Applied Geodynamics: A Numerical
Study.&rdquo; *Geoch. Geoph. Geosystems* 12: Q07004/1&ndash;26.

</div>

</div>

  [benchmarks/solkz/]: benchmarks/solkz/
  [2]: #fig:solkz
  [3]: #fig:solkz2
