# Kinematically-driven 2d oceanic subduction

*This section was contributed by Anne Glerum and Yingying Li.*

This subduction example is based on the benchmark effort of Quinquis et al.,
of which initial results were published in {cite:t}`quinquis:2014`. In four
increasingly complex cases we will go from isoviscous materials without any
temperature effects to a fully thermo-mechanical, nonlinear, strain-weakened
visco-plastic, externally-driven model of oceanic subduction. The setup is
outlined in {numref}`fig:QQ_setup`. The models are run for 15 My
and slab tip depth, trench location, root-mean-square (RMS) velocity and temperature, and viscous
dissipation are monitored. In addition, we will discuss the effects of the element
size of the subduction interface and crustal layers, viscosity averaging and
the solver tolerance.

```{figure-md} fig:QQ_setup
<img src="setup_Quinquis2014.*" style="width:80.0%" />

 Case 4 model setup. Copied from {cite:t}`quinquis:2014`.
```

## Case 1: Simple rheology

The Case 1 model setup considers seven materials (compositional fields) apart
from the background sublithospheric mantle (see {numref}`fig:QQ_case1_setup`), including (1) the Bulk Oceanic Composition (BOC), (2) Serpentinized HarzBurgite
(SHB) and (3) "thermal" layer of the overriding plate (OP), (4) the BOC, (5) SHB and (6)
"thermal" layer of the subducting plate (SP), and (7) the weak seed.

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
velocity of the slab and the whole model domain, and the viscous
dissipation in the slab and total model domain. The computation of trench position
is through a plugin called 'trench_location'. It needs to be compiled as a shared library
that can be called from the input file.

```{literalinclude} Case1_postprocessing.prm

```


```{figure-md} fig:QQ_case1_setup
<img src="Case1_t0.*" style="width:60.0%" />

 Case 1 density, viscosity and velocity at time zero.
```

We run the Case 1 setup for 15 My of model time. The diagnostic quantities in
{numref}`fig:QQ_case1_diagnostics` show three stages of model evolution: first trench advance (top right plot), then free subduction (increasing slab RMS velocity), and
after about 13 My interaction between the slab and bottom boundary at 660 km
depth, which slows down the slab. The slab then curves inward along the bottom
boundary. This can also be seen in {numref}`fig:QQ_case1_results`.

```{figure-md} fig:QQ_case1_diagnostics
<img src="Case1_diagnostics.*" />

 Case 1 diagnostic quantities of ASPECT
```


```{figure-md} fig:QQ_case1_results
<img src="Case1_visc.*" style="width:60.0%" />

 Case 1 viscosity snapshots at 8 and 15 My.
```
## Case 2a: Model with thermal evolution

Case 2a builds on Case 1 by including a heterogeneous initial temperature field
that combines the plate cooling model for the temperature in the plates with an
adiabatic temperature gradient below the plates. For the plate cooling model,
we assume that the overriding plate is 40 My old and the subducting plate 70 My.
The temperature function (including a summation term) is too complex to write out as a function
expression in the input file, so it is put in a new initial temperature plugin
called ’subduction plate cooling’. This plugin needs to be compiled as a shared
library that can be called from the input file:

```{literalinclude} Case2a_temperatures.prm

```
The initial temperature distribution can be seen in  {numref}`fig:QQ_case2a_temperature_t0`

```{figure-md} fig:QQ_case2a_temperature_t0
<img src="Case2a_temp_t0.*" />

 Case 2a temperature distribution at time 0
```

The thermal expansivity is set to zero in Case 2a, so that the thermal evolution is
not coupled to the mechanical evolution through the density. A high conductivity of
k = 183.33 Wm<sup>-1</sup>K<sup>-1</sup> is assigned to the upper mantle to maintain the mantle adiabat.

```{literalinclude} Case2a_materialmodel.prm

```

We monitor the temperature evolution by tracking the
root-mean-square (RMS) temperature and the isotherm depth using
two plugins called "composition_trms_statistics" and "isotherm depth".
To this end, we add the postprocessors "composition RMS temperature statistics"
and "isotherm depth" to the prm file. The RMS temperature
of the subducting plate (comprising the BOC_SP, SHB_SP, thermal_SP)
as well as the whole domain is tracked over time, as is the
800 degrees Celsius isotherm.

```{literalinclude} Case2a_postprocessing.prm

```

{numref}`fig:QQ_case2a_diagnostics` shows the RMS temperature,
800 degrees Celsius isotherm depth as well as the other diagnostics. For reference,
results from the ELEFANT software {cite}`thieulot:2015` are also shown. The 800 degrees Celsius isotherm
depth continuously increases over time, but towards the end at a slower rate.
The RMS temperature of the whole domain keeps decreasing due to the
cold subducting plate entering the model domain. The slab RMS temperature
decreases at the initial stage, then it goes up due to the slab being warmed up by the
mantle. Note that the other diagnostics are indeed the same as for Case 1.

```{figure-md} fig:QQ_case2a_diagnostics
<img src="Case2a_diagnostics.*" />

 ASPECT Case 1 and Case 2a as well as Elefant Case 2a diagnostic quantities
```

## Case 2b: Model with a temperature-dependent density

Case 2b builds on Case 2a.  Compared to Case 2a, the thermal expansivity in Case 2b
is set to 2.5e-5 K<sup>-1</sup> instead of 0 K<sup>-1</sup>. This renders the density
(and therefore the Stokes equations) temperature dependent. In addition, the benchmark uses different reference temperatures for each composition, but only one value can be specified for the reference temperature
in ASPECT. Therefore, we adapted the reference densities for the different compositional fields to reflect the different reference temperatures. All the other parameters are
identical to Case 2a:

```{literalinclude} Case2b_materialmodel.prm

```
The density evolution can be seen in {numref}`fig:QQ_case2b_density_evolution`. As the slab is heated up by the surrounding mantle, its density contrast with the mantle decreases, and so does its negative buoyancy.

```{figure-md} fig:QQ_case2b_density_evolution
<img src="Case2b_density_evolution.*" />

 Case2b density evolution
```

Therefore, subduction occurs at a slower pace than in Case 1 and 2a (see {numref}`fig:QQ_case2b_diagnostics`).
The subducting plate does not reach the bottom in 15 Myr, i.e., the third stage of Case 1 and 2a is missing.

```{figure-md} fig:QQ_case2b_diagnostics
<img src="Case2b_diagnostics.*" />

 ASPECT Case2b diagnostic quantities
```

The total RMS velocity is however higher than in the previous cases (see {numref}`fig:QQ_case2b_diagnostics`) as an additional convection cell to the left of the slab increases the sub-lithospheric velocities (see {numref}`fig:QQ_case2b_velocity_evolution`).

```{figure-md} fig:QQ_case2b_velocity_evolution
<img src="Case2b_velocity_evolution.*" />

 Case2b velocity evolution
```

Case 2b remains relatively simple, as all materials are linear viscous. In Case 3 (to be added) we will address nonlinear viscoplastic rheologies.
