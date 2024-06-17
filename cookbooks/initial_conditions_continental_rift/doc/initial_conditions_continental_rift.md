# Initial and boundary conditions for continental rifting

*This section was contributed by Anne Glerum.*

Continental rifting is a process that thins the continental
lithosphere and can eventually lead to full break-up of that lithosphere
and seafloor spreading at a mid-ocean ridge {cite}[e.g., ][]`brune:2023`.
We can simulate continental rifting under far-field extensional
forces by prescribing certain velocities or stresses on the lateral
boundaries of our model domain (as in, for example, section 
{ref}`sec:cookbooks:continental_extension`).
Here we choose to prescribe outward horizontal velocities on
the lateral boundaries.
In order to localize deformation under this prescribed 
extension (and avoid localizing at the model boundaries), 
we have to introduce initial weaknesses to the
model setup. The plugins prescribed in this cookbook set up
a lithospheric (upper crust, lower crust and mantle lithosphere) and asthenosphere
composition and temperature structure and then
add weaknesses by 1) smoothly perturbing the lithosphere
at user-specified locations and by 2) adding initial plastic
strain around these same locations. This strain represents
inherited weakness from for example prior deformation phases.
To represent other deviations of the reference lithosphere,
such as cratonic lithosphere, 
additional regions of different lithospheric layer thicknesses
can be specified.
In case the perturbations of the lithosphere in combination with
a free surface or another type of mesh deformation induce large initial 
velocities, an initial topography plugin approximates the topography resulting from
isostatic equilibrium with respect to the reference (unperturbed) lithosphere.

Apart from prescribing tensional velocities on the lateral boundaries of
the box domain, we use the surface processes software FastScape to
modify the top boundary and a traction boundary condition for the bottom
boundary. This allows for the inclusion of sediment erosion, transport
and deposition, and the automatic balance of mass within the model domain.

## Initial lithosphere composition and temperature

```{figure-md} fig:setup
<img src="Initial_density_with_craton.*" style="width:80.0%" />

Lithospheric structure with perturbations that represent an
incipient rift and a cratonic region.
```
To set up the reference lithospheric structure, we specify the Layer thicknesses
of the upper crust, lower crust and mantle lithosphere within the Lithosphere with rift 
subsection (see input-file snippet below). The perturbations of these lithospheric layers - representing an
incipient rift and localizing deformation - follow Gaussian curves. 
We can either thicken or thin the layers, depending on
the sign of each Gaussian amplitude specified: a negative amplitude thickens
the corresponding layer and vice versa. The width of the perturbation
is determined by the sigma of the Gaussians. The centre of the Gaussian
curve will be the rift centre and is set with the parameter Rift axis line segments: 


```{literalinclude} initial_composition_lithosphere.prm
```

Note that you can specify multiple segments, both in 2D and in 3D. By changing
the thickness of the layers, we simultaneously change the temperature profile,
and thus the strength distribution, creating weaker lithosphere.

In addition to an incipient rift, other deviations from the reference
lithospheric thicknesses can also be prescribed (e.g., representing thick, cold cratonic
lithosphere). In the above snippet,
these are specified as polygons (multiple polygons are allowed) with
their own layer thicknesses. The transition from reference lithosphere
to polygon lithosphere is implemented as a hyperbolic tangent with a
centre point defined by the coordinates of the polygon and a half-width (10 km in this case).


```{literalinclude} initial_temperature.prm
```


For the initial temperature we combine two plugins: one that prescribes
a continental steady-state geotherm (lithosphere with rift) and one that computes an adiabatic
temperature profile (see snippet above). The geotherm implementation is based on Chapter 4.6
of Turcotte and Schubert. For each x-coordinate, it solves a steady-state
heat conduction equation taking into account the respective densities,
heat productivities and thermal conductivities of the three lithospheric
layers and the user-set temperature at the top (Surface temperature) and 
bottom of the lithosphere (LAB isotherm temperature). Below the lithosphere a NaN is returned, such that we can replace
this NaN by an adiabatic temperature profile. Care has to be taken to match
the user-specified LAB temperature and the adiabatic temperature at the same
depth by modifying the adiabatic surface temperature. One can choose to prescribe
the LAB temperature up to a certain compensation depth (e.g., the bottom of the
craton), such that the LAB and adiabatic temperatures match everywhere (Fig. {numref}`fig:initial_temperature`).

```{figure-md} fig:initial_temperature
<img src="Initial_temperature_with_craton.*" style="width:60.0%" />

 The initial temperature distribution.
```
## Inherited heterogeneity

```{literalinclude} initial_composition_strain.prm
```

To include inherited heterogeneity we also include initial plastic strain
around the rift centre (Fig. {numref}`fig:initial_strain`). This requires a
compositional field called plastic_strain and the plugin rift box initial plastic strain.
We generate random strain values between 0 and 1 that
are subsequently multiplied with a Gaussian distribution for which we specify
the maximum amplitude and sigma. To control the depth extent of the initial
strain, we specify a depth around which the amplitude is smoothed out to zero
with a hyperbolic tangent. The initial strain thus produced will be the same
for each run with the same specifications if the same number of MPI processes is used.
The plastic strain can be used to weaken the cohesion and/or internal angle of friction
of the plastic rheology, thus again providing weaknesses that localize deformation.

```{figure-md} fig:initial_strain
<img src="Initial_plastic_strain_with_craton.*" style="width:60.0%" />

 The initial plastic strain distribution.
```

```{literalinclude} initial_topography.prm
```

## Initial topography

Sometimes the combination of a free surface with lateral variations in the
lithospheric thicknesses leads to very high velocities in the first timesteps.
To avoid a drunken sailor effect and/or very small timesteps, we can then
prescribe initial topography that roughly fulfills isostatic equilibrium (Fig. {numref}`fig:initial_topography`).
The required topography is computed using the reference density for each material, 
so it neglects temperature effects on the density.


```{figure-md} fig:initial_topography
<img src="Initial_topography_with_craton.*" style="width:60.0%" />

 The initial topography that roughly isostatically balances the rift perturbation and craton region
 with respect to the reference lithosphere.
```

## Boundary conditions

```{literalinclude} velocity_boundary_conditions.prm
```
Extension is driven by outward horizontal velocities on the left and
right boundary of the model domain. To compensate for this outflow,
we open the bottom boundary to flow by prescribing a traction. The 
prescribed traction is computed as the initial lithostatic pressure 
at a user-specified point. We pick this point to lie in the reference
(unperturbed) lithosphere close to the left domain boundary.

For the top boundary, we pick the mesh deformation plugin fastscape,
which provides an interface to the FastScape landscape evolution model.
This plugin computes changes to the mesh's surface topography based on
river incision, hill slope diffusion, sediment deposition and marine diffusion. 
For all the different parameters relating to FastScape, have a look at 
section {ref}`sec:cookbooks:fastscape_eroding_box`. We allow the top corner
of the right boundary to move with the surface, but not the left top corner point,
as this is our reference location.
