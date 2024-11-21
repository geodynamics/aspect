(sec:cookbooks:mantle-conv-annulus)=
# Mantle convection with continents in an annulus

*This section was contributed by Cedric Thieulot and Erik van der Wiel.*

**UNDER CONSTRUCTION: This cookbook does currently not produce the results shown
in the figures. The reason for that is that the original model of
{cite:t}`vanderwiel:etal:2024a` used a more complicated parameter file and
a different version of ASPECT. Do not rely on the results of this cookbook
at the moment.**

The setup for this experiment originates in {cite:t}`vanderwiel:etal:2024a`.
The domain is an annulus with an inner radius of 3480 km and an outer radius of 6370 km.
The Isentropic Compression Approximation (ICA) is used {cite}`gassmoller:etal:2020`, which is
the default approximation for compressible flow in ASPECT.

Flow in the model is governed by the visco-plastic flow equations
that describe dislocation and diffusion creep and we use the Drucker-Prager
yield criterion to limit viscous stresses {cite}`glerum:etal:2018`.
The grain size in the diffusion creep is assumed constant and therefore
incorporated in the pre-factor A.
Note that in the lower mantle flow is based on diffusion creep only.

Both the inner and outer boundaries are free-slip boundaries so
there is no external kinematic forcing on the model. The resulting
existing rotational null space is removed by imposing no-net-rotation of
the mantle. The boundaries have fixed temperatures of 3700 K and 300 K.
Within the domain two compositional fields are defined: mantle and
continent. The mantle domain consists of three different regions separated
by the two major phase changes occurring at ~410- and ~660 km
depths. We use the geodynamic World Builder {cite}`Fraters2019c` to
set the initial temperature and compositional field distribution
(see {numref}`fig:init_visc`, {numref}`fig:init_compositions` and  {numref}`fig:init_temperature`).


```{figure-md} fig:init_visc
<img src="viscosity.*" width="90%" />

Initial viscosity field.
```

```{figure-md} fig:init_compositions
<img src="continents.*" width="90%" />

The initial position of the three continents.
```

```{figure-md} fig:init_temperature
<img src="temperature.*" width="90%" />

Initial temperature field.
```

Similar to {cite:t}`ulvrova:etal:2019`, we use viscous rafts as "continents" to aid with
modelling one-sided subduction systems. We use three continents of
5000, 5000 and 3000 km length, covering roughly 30% of the model
surface. To avoid subduction of the continents, we assign the viscous
rafts a reference density of 2916 kg/m<sup>3</sup>, which corresponds to a 400
kg/m<sup>3</sup> density contrast with the reference mantle density.
Furthermore, to avoid deformation of the continents, we assign
them a viscosity that equals the maximum cut-off viscosity in the model,
i.e. 6e24 Pa s. We also configure three subduction zones
in the initial set-up to seed the model with initial negative
buoyancy.

As the goal of the cookbook (and the corresponding publication in {cite:t}`vanderwiel:etal:2024a`) is to investigate
how average slab sinking rates relate to the vigor of mantle convection and mixing,
various rheologies are considered in the lower mantle, as
shown in {numref}`fig:visc-profiles`.
We will focus on the 'R' (reference) profile in this cookbook (green line).

```{figure-md} fig:visc-profiles
<img src="viscosity_profiles.*" width="90%" />

Radial viscosity profiles for the three models discussed in the original
publication (solid lines)
compared to viscosity profiles as used in other studies.
See original publication for the references corresponding to the
dashed lines.
```

After 500 Myr we see multiple subductions taking place as well as
many plumes/upwellings, as shown in {numref}`fig:evolution`

```{figure-md} fig:evolution
<img src="evolution.*" width="90%" />

Results of the reference model (R) after 500 Ma of mantle convection simulation. Snapshots of the model at t = 500 Ma are shown for the viscosity (a) and
temperature (b) within the modelled domain, showing five slabs actively subducting below the continents (pink).
(c) Solid lines indicate average surface velocity (blue)
and mobility M (red) of model R as well as their dashed time-averages Vsurf=2.1 cm/a and M = 2.218.
Mobility is defined as the ratio of rms surface velocity to rms velocity averaged over the entire 3D domain {cite}`tackley:2000`.
(d) Radial velocity of all tracers defined as
slabs plotted at their position in the model where blue indicates sinking slabs and red slowly rising. Shown tracers are 300 K colder than the radial averaged
temperature at similar depth and automatically obtained (see methods section of the publication).
(e) 2D-histogram showing the tracer depth vs the radial velocity in 25 km by 2 mm/a bins. The
colour indicates the number of particles within a certain bin.
```
