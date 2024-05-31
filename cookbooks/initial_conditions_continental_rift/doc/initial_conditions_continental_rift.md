# Initial conditions for continental rifting

*This section was contributed by Anne Glerum.*

In order to localize deformation under prescribed far-field
extension, initial weaknesses have to be introduced to the
model setup. Here we discuss plugins that 1) perturb the lithosphere
consisting of upper crust, lower crust and mantle lithosphere
at a user-specified location and 2) add initial plastic
strain is around this same location. This strain represents
inherited weakness for example from prior deformation phases.
To represent other deviations of the reference lithosphere,
additional regions of different lithospheric layer thicknesses
can also be specified.

The initial temperature plugin prescribes a steady-state
continental geotherm within the lithosphere and an adiabat
below. The temperature field is perturbed in the same way
as the lithosphere compositional structure.

In case the perturbations of the lithosphere in combination with
the free surface induce large velocities initially, an initial
topography plugin approximates the topography resulting from
isostatic equilibrium.


## Initial lithosphere compositional fields and temperature

```{figure-md} fig:setup
<img src="Initial_density_with_craton.*" style="width:80.0%" />

Lithospheric structure with perturbations that represent an
incipient rift and a cratonic region.
```

The perturbation of the lithospheric layers that represents an
incipient rift follows a Gaussian
curve. We can either thicken or thin the layers, depending on
the sign of the Gaussian amplitude: a negative amplitude thickens
the corresponding layer and vice versa. The width of the perturbation
is determined by the sigma of the Gaussian. The centre of the Gaussian
curve will be the rift centre. The parameters that need to be set are:


```{literalinclude} initial_composition_lithosphere.prm
```

In addition to a rift centre, other deviations from the reference
lithospheric thicknesses can also be prescribed. In the above snippet,
these are specified as polygons (multiple polygons are allowed) with
their own layer thicknesses. The transition from reference lithosphere
to polygon lithosphere is implemented as a hyperbolic tangent with a
centre point and a half-width.


```{literalinclude} initial_temperature.prm
```

For the initial temperature we combine two plugins: one that prescribes
a continental steady-state geotherm and one that computes an adiabatic
temperature profile. The geotherm implementation is based on Chapter 4.6
of Turcotte and Schubert. For each x-coordinate, it solves a steady-state
heat conduction equation taking into account the respective densities,
heat productivities and thermal conductivities of the three lithospheric
layers. Below the lithospheric a NaN is returned, such that we can replace
this NaN by an adiabatic temperature profile. Care has to be taken that the
adiabatic surface temperature matches the LAB temperature.

```{figure-md} fig:initial_temperature
<img src="Initial_temperature_with_craton.*" style="width:60.0%" />

 The initial temperature distribution.
```

```{literalinclude} initial_composition_strain.prm
```

To include inherited heterogeneity we also include initial plastic strain
in the rift centre (Fig. {numref}`fig:initial_strain`).
We generate random strain values between 0 and 1 that
are subsequently multiplied with a Gaussian distribution for which we specify
the maximum amplitude and sigma. To control the depth extent of the initial
strain, we specify a depth around which the amplitude is smoothed to zero
with a hyperbolic tangent.

```{figure-md} fig:initial_strain
<img src="Initial_plastic_strain_with_craton.*" style="width:60.0%" />

 The initial plastic strain distribution.
```

```{literalinclude} initial_topography.prm
```

Sometimes the combination of a free surface with lateral variations in the
lithospheric thicknesses leads to very high velocities in the first timesteps.
To avoid a drunken sailor effect and/or very small timesteps, we can then
prescribe initial topography that roughly fulfills isostatic equilibrium (Fig. {numref}`fig:initial_topography`).
The required topography is computed using the reference density, so it neglects
temperature effects on the density.


```{figure-md} fig:initial_topography
<img src="Initial_topography_with_craton.*" style="width:60.0%" />

 The initial topography that roughly isostatically balances the rift perturbation and craton region
 with respect to the reference lithosphere.
```
