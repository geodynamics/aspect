#### Kinematically-driven 2d oceanic subduction

*This section was contributed by Anne Glerum.*

This subduction example is based on the benchmark effort of Quinquis et al.,
of which initial results were published in (Quinquis 2014). In four
increasingly complex setups we will go from isoviscous materials without any
temperature effects to a fully thermo-mechanical, nonlinear, strain-weakened
visco-plastic, externally-driven model of oceanic subduction. The setup of the
most complex case is outlined in Fig.&nbsp;[1][]. The models are run for 15 My
and slab tip depth, trench location, RMS velocity and temperature, and viscous
dissipation are monitored. In addition, we discuss the effects of the element
size of the subduction interface and crustal layers, viscosity averaging and
the solver tolerance.

<figure>
<img src="cookbooks/kinematically_driven_subduction_2d/doc/setup_Quinquis2014.png" id="fig:QQ_setup" style="width:80.0%" alt="Case 4 model setup. Copied from Quinquis (2014)." /><figcaption aria-hidden="true">Case 4 model setup. Copied from Quinquis (2014).</figcaption>
</figure>

##### Case 1: Simple geometry and rheology

The Case 1 model setup considers three materials (compositional fields) apart
from the background sublithospheric mantle (see Fig.&nbsp;[4][]):

1.  the lithosphere of the overriding plate (combining the BOC, SHB and
    thermal layer of the overriding plate in Fig.&nbsp;[1][]),

2.  the crust of the subducting plate (weak seed and BOC combined), and

3.  the mantle lithosphere of the subducting plate (SHB and thermal layer
    combined).

The geometry of these compositions is implemented as follows:

``` prmfile
```

No differences in material properties exist, except for density and viscosity,
so we use the multicomponent material model. To keep the boundaries between
fields as sharp as possible in terms of viscosity, we use the maximum
composition to determine the viscosity in each evaluation point (note that
this can be harder on the solver):

``` prmfile
```

Temperature effects are ignored. Subduction is driven by prescribed in- and
outflow through the right boundary (with a gradual transition of the flow
direction), all other boundaries are free slip. The volume of material that
flows in is balanced by the volume that is prescribed to flow out (this is
important as the model is incompressible). A weak crust along the plate
interface and the subducting lithosphere facilitates subduction. Through the
function plugin, we prescribe the in- and outflow:

``` prmfile
```

To follow the slab on its descent into the mantle, we use adaptive mesh
refinement based on viscosity in combination with the minimum refinement
strategy to ensure sufficient resolution in the crust and weak zone that allow
the slab to detach from the surface:

``` prmfile
```

To monitor the model evolution, several diagnostic quantities are tracked over
time: the depth of the tip of the slab, the position of the trench, the RMS
velocity of the slab and the whole model domain, and the work done (viscous
dissipation) in the slab and total model domain. The computation of these
quantities in done in several new plugins:

``` prmfile
```

<img src="cookbooks/kinematically_driven_subduction_2d/doc/Case1_dens_t0.png" title="fig:" id="fig:QQ_case1_setup" style="width:60.0%" alt="Case 1 density, viscosity and velocity at time zero." />
<img src="cookbooks/kinematically_driven_subduction_2d/doc/Case1_visc_t0.png" title="fig:" id="fig:QQ_case1_setup" style="width:60.0%" alt="Case 1 density, viscosity and velocity at time zero." />
<img src="cookbooks/kinematically_driven_subduction_2d/doc/Case1_vel_t0.png" title="fig:" id="fig:QQ_case1_setup" style="width:60.0%" alt="Case 1 density, viscosity and velocity at time zero." />

We run the Case 1 setup for 15 My of model time. The diagnostic quantities in
Fig.&nbsp;[5][] show three stages of model evolution: first trench advance
(top right plot), then free subduction (increasing slab RMS velocity), and
after about 13 My interaction between the slab and bottom boundary at 660 km
depth, which slows down the slab. The slab then curves inward along the bottom
boundary. This can also be seen in Fig&nbsp;[7][].

<figure>
<img src="cookbooks/kinematically_driven_subduction_2d/doc/Case1_diagnostics.*" id="fig:QQ_case1_diagnostics" alt="Case 1 diagnostic quantities of ASPECT, Sulec and Elefant results." /><figcaption aria-hidden="true"><em>Case 1 diagnostic quantities of ASPECT, Sulec and Elefant results.</em></figcaption>
</figure>

<img src="cookbooks/kinematically_driven_subduction_2d/doc/Case1_visc_8My.png" title="fig:" id="fig:QQ_case1_results" style="width:60.0%" alt="Case 1 viscosity snapshots at 8 and 15 My." />
<img src="cookbooks/kinematically_driven_subduction_2d/doc/Case1_visc_15My.png" title="fig:" id="fig:QQ_case1_results" style="width:60.0%" alt="Case 1 viscosity snapshots at 8 and 15 My." />

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Quinquis2014" class="csl-entry">

Quinquis, M. 2014. &ldquo;<span class="nocase">A numerical study of subduction
zone dynamics using linear viscous to thermo-mechanical model setups including
(de)hydration processes</span>.&rdquo; PhD thesis, Charles University.

</div>

</div>

  [1]: #fig:QQ_setup
  [4]: #fig:QQ_case1_setup
  [5]: #fig:QQ_case1_diagnostics
  [7]: #fig:QQ_case1_results
