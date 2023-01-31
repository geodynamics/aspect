# Kinematically-driven 2d oceanic subduction

*This section was contributed by Anne Glerum.*

This subduction example is based on the benchmark effort of Quinquis et al.,
of which initial results were published in {cite:t}`quinquis:2014`. In four
increasingly complex setups we will go from isoviscous materials without any
temperature effects to a fully thermo-mechanical, nonlinear, strain-weakened
visco-plastic, externally-driven model of oceanic subduction. The setup of the
most complex case is outlined in {numref}`fig:QQ_setup`. The models are run for 15 My
and slab tip depth, trench location, RMS velocity and temperature, and viscous
dissipation are monitored. In addition, we discuss the effects of the element
size of the subduction interface and crustal layers, viscosity averaging and
the solver tolerance.

```{figure-md} fig:QQ_setup
<img src="setup_Quinquis2014.*" style="width:80.0%" />

 Case 4 model setup. Copied from {cite:t}`quinquis:2014`.
```

## Case 1: Simple geometry and rheology

The Case 1 model setup considers three materials (compositional fields) apart
from the background sublithospheric mantle (see {numref}`fig:QQ_case1_setup`):

1.  the lithosphere of the overriding plate (combining the BOC, SHB and
    thermal layer of the overriding plate in {numref}`fig:QQ_setup`),

2.  the crust of the subducting plate (weak seed and BOC combined), and

3.  the mantle lithosphere of the subducting plate (SHB and thermal layer
    combined).

The geometry of these compositions is implemented as follows:

```{literalinclude} Case1_compositions.prm
```

No differences in material properties exist, except for density and viscosity,
so we use the multicomponent material model. To keep the boundaries between
fields as sharp as possible in terms of viscosity, we use the maximum
composition to determine the viscosity in each evaluation point (note that
this can be harder on the solver):

```{literalinclude} Case1_materialmodel.prm
```

Temperature effects are ignored. Subduction is driven by prescribed in- and
outflow through the right boundary (with a gradual transition of the flow
direction), all other boundaries are free slip. The volume of material that
flows in is balanced by the volume that is prescribed to flow out (this is
important as the model is incompressible). A weak crust along the plate
interface and the subducting lithosphere facilitates subduction. Through the
function plugin, we prescribe the in- and outflow:

```{literalinclude} Case1_velocity.prm
```

To follow the slab on its descent into the mantle, we use adaptive mesh
refinement based on viscosity in combination with the minimum refinement
strategy to ensure sufficient resolution in the crust and weak zone that allow
the slab to detach from the surface:

```{literalinclude} Case1_meshrefinement.prm
```

To monitor the model evolution, several diagnostic quantities are tracked over
time: the depth of the tip of the slab, the position of the trench, the RMS
velocity of the slab and the whole model domain, and the work done (viscous
dissipation) in the slab and total model domain. The computation of these
quantities in done in several new plugins:

```{literalinclude} Case1_postprocessing.prm
```


```{figure-md} fig:QQ_case1_setup
<img src="Case1_t0.*" style="width:60.0%" />

 Case 1 density, viscosity and velocity at time zero.
```

We run the Case 1 setup for 15 My of model time. The diagnostic quantities in
{numref}`fig:QQ_case1_diagnostics` show three stages of model evolution: first trench advance
(top right plot), then free subduction (increasing slab RMS velocity), and
after about 13 My interaction between the slab and bottom boundary at 660 km
depth, which slows down the slab. The slab then curves inward along the bottom
boundary. This can also be seen in {numref}`fig:QQ_case1_results`.

```{figure-md} fig:QQ_case1_diagnostics
<img src="Case1_diagnostics.*" />

 Case 1 diagnostic quantities of ASPECT, Sulec and Elefant results.
```


```{figure-md} fig:QQ_case1_results
<img src="Case1_visc.*" style="width:60.0%" />

 Case 1 viscosity snapshots at 8 and 15 My.
```
