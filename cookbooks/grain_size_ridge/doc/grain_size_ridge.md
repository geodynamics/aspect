(sec:cookbooks:grain-size-ridge)=
# Sublithospheric convection beneath an oceanic plate with a grain size dependent rheology

*This section was contributed by Juliane Dannberg.*

The input file for this model can be found at
[cookbooks/grain_size_ridge/grain_size_ridge.prm](https://www.github.com/geodynamics/aspect/blob/main/cookbooks/grain_size_ridge/grain_size_ridge.prm)

Note that this model may take a significant time to run due to the high
resolution and large number of time steps. We recommend to either run it
on a workstation or cluster or to reduce the model resolution.

This model features the spreading of an oceanic plate away from a mid-ocean
ridge, with the ridge being in the top left corner of the model. Its vertical
extent is 410 km and and its horizontal extent is 4000 km. Since the plate
moves with a velocity of 5 cm/yr, this leads to a maximum plate age of 80 Myr.
Material flows in from the bottom, and then leaves the domain at the right side.
The initial temperature is the adiabatic temperature of the mantle with an
added cold boundary layer becoming thicker with increasing distance from the
ridge (increasing plate age). The mineral grain size starts from a constant
value of 5 mm and then evolves over time. Grains grow faster at higher
temperature, and the grain size is reduced by deformation in the dislocation
creep regime.


```{figure-md} fig:grain-size-ridge
<img src="setup.png" style="width:60.0%" />

 Setup of cooling oceanic plate model with grain size dependent rheology. Top: Background colors show temperature, streamlines illustrate the flow. Bottom: viscosity.
```

## Grain size evolution.

To model the evolution of the mineral grain size and to take this grain size
into account in our rheological model, we use the grain size material model,
where these equations are implemented. This model implements a composite
diffusion--dislocation rheology, with diffusion creep depending on the
evolving grain size, and deformation by dislocation creep reducing the grain
size. For more details on this material model and how to use it for global
mantle convection models, see {cite:t}`dannberg:grain:size`.
We have to specify the parameters for both creep mechanisms and for the grain
size evolution.

```{literalinclude} material_model.part.prm
```

Many of these parameters, in particular the prefactors and activation volumes,
have large uncertainties. The values used here should therefore only be seen
as example values. It often makes sense to choose their values based on the
viscosity profile one wants to use. These viscosity profiles are usually
constrained by observations, such as the geoid or mineral physics data, and
therefore provide additional constraints to the values we get from experiments.

Note that diffusion creep depends on the third power of the grain size as
specified by the `Diffusion creep grain size exponent` parameter, and that
the stress-dependence of dislocation creep is controlled by the
`set Dislocation creep exponent` parameter.


## Tracking the grain size on particles.

To track this evolving mineral grain size on particles, we have to have a
compositional field with the name `grain_size`, we have to specify that
we want to advect this field using the particle method, and we have to select
`grain_size` as the particle property to be mapped to this compositional field:

```{literalinclude} fields.part.prm
```

In addition, we have to specify a number of different parameters controlling
how to initialize and advect these particles, how to interpolate the computed
values back to the mesh and how to control the load balancing.
We here choose to add and remove particles in such a way that there are always
at least 40 and not more than 640 particles in each cell. This means that a
cell with 160 particles can either be refined once (resulting in about 40
particles in each new cell) or coarsened once (resulting in about 640 in the
new bigger cell) without having to add or remove particles. We also use the
bilinear least squares interpolation scheme to interpolate from the particles
to the grid, which is more accurate than schemes based on, for example, the
cell average. To make sure this scheme does not generate any over- or under-
shoots in the grain size, we switch on the linear least squares limiter. To
use this scheme, we also have to update the ghost particles, because the
interpolation scheme requires knowledge about particles in neighboring cells.
Finally, we use a second Order Runge Kutta integrator to advect the particles.

```{literalinclude} particles.part.prm
```

## Model evolution.

We run the model for 200 million years. In the first few tens of millions of years
there are no large changes in the temperature structure, except that the
thermal boundary layer increases slightly in thickness. In addition, the grain
size evolves from its initially constant value, with the grain size being
reduced due to deformation. This leads to a band of particularly low grain
size in the lowermost part of the oceanic plate. Grain size is also reduced
in the asthenesphere, but to a lesser degree.

At 58.8 Myr, the first downwelling from the base of the lithosphere starts to
form. Due to the strong temperature-dependence of the viscosity, only the
lowermost part of the plate participates in this downwelling, which is
therefore very thin compared to the plate. Because of the strong deformation
within the downwelling (which reduces the grain size) and its lower
temperature (which hinders grain growth), the grain size in these convective
instabilities is smaller compared to the surrounding mantle. Over time,
the downwelling moves towards the left and younger plate ages, and additional
downwellings start to develop until the model reaches a quasi-steady state.

```{figure-md} fig:grain-size-ridge-result
<img src="evolution.svg" style="width:60.0%" />

Evolution of the temperature (top panels) and grain size (bottom panels) in the model.
The first downwelling starts at 58.8 Myr. The full model evolution can be found
on [YouTube](https://youtu.be/qeeucnS7A58).
```

Note that for simplicity, our model uses a constant thermal conductivity of
4.7 W/m/K, whereas the thermal conductivity is likely lower in the lithosphere
and is also temperature-dependent. Therefore, our lithosphere is thicker (and
downwellings occur earlier) than on Earth.
