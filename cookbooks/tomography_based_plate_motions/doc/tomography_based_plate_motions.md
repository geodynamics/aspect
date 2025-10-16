```{tags}
category:cookbook
feature:2d
feature:spherical
feature:compositional-fields
feature:data-integration
feature:compressibility
feature:nonlinear-solver
```

(sec:cookbooks:tomography-plate-motions)=
# Mantle convection using tomography data

*This section was contributed by Arushi Saxena, Juliane Dannberg and Ren&eacute; Gassm&ouml;ller.*

Mantle flow models based on present-day observational constraints can help us to investigate several global
and regional tectonic processes. This cookbook describes how to set up such a model based on input tomography
data and plate boundary locations. These models are typically computationally very expensive in 3D to make
meaningful comparisons with the observations, therefore, we explain the workflow for a 2D spherical shell.

## Description of input files

In this cookbook, the input tomography data uses the LLNL-G3D-JPS model {cite}`simmons:etal:2019`, available
at the [IRIS DMC webpage](http://ds.iris.edu/ds/products/emc-earthmodels/), read in as an ASCII compositional data file.
The input file is located in `cookbooks/tomography_based_plate_motions/input_data/` and has the following format
```{literalinclude} input_tomography.part.txt
```
where grain_size is set to a constant value of 5mm, Vp and Vs are the P- and S-wave seismic velocities, respectively,
Vs_anomaly is the S-wave velocity anomaly relative to a reference velocity profile, faults represent the composition
value at plate boundaries, and cratons represent the composition value within the continental cratons. Both faults and
cratons have composition values set to 0 in this ASCII file because they are defined separately using the Geodynamics
World Builder (GWB) {cite}`Fraters2019c` and the input file `cookbooks/tomography_based_plate_motions/input_file/world_builder_smac_cratons_faults_2D.json`.

For defining plate boundary geometry, we use the Nuvel model {cite}`demets:etal:1990` input as "fault" features and
the locations of cratons taken from {cite}`nataf:etal:1996` input as "continental plate" features in GWB {cite}`Fraters2019c`.
The compositional value of faults transitions smoothly from 1 at the fault trace to 0 over a width of 50 km on either
side of the fault. We do this to ensure that the prescribed plate boundary viscosity which is based on the fault
composition is also smooth. Note that while the input file includes plate boundary geometry and cratons for the whole Earth,
we only compute the compositions along a cross-section in this cookbook. The input faults and cratons in a 3D spherical model
would look like {numref}`fig:composition`.

```{figure-md} fig:composition
<img src="Fig_composition.*" />

The composition of plate boundary geometry ("faults") and cratons ("continents") using the input file coordinates for a 3D model.
```

## Computed material properties
A common practice in mantle convection models using tomography data is to use a different model for the uppermost mantle
structure because of the limited seismic resolution to accurately resolve the lithospheric structure and slabs. For this cookbook,
we use a temperature model, TM1 {cite}`osei:etal:2018`, located in `cookbooks/tomography_based_plate_motions/input_data/upper_mantle_TM1_2D.txt`.
The user can control the depth until which we want to use this temperature model using the `uppermost mantle thickness` parameter.
{numref}`fig:Vs_temperatures` shows the input Vs anomalies and the computed temperatures in the model.

```{figure-md} fig:Vs_temperatures
<img src="Fig_tomography_temperature.*" />

The input Vs anomalies (left) and the computed temperatures (right). Note that we only use tomography-based temperatures below the
uppermost mantle thickness depth, marked by the white contour line.
```

We compute the density from the thermal anomalies in the TM1 model relative to a reference temperature of 293 K as done by
{cite:t}`osei:etal:2018`. We use constant thermal expansion coefficients and compressibilities within the crust, lithospheric mantle and
asthenosphere, respectively. To define the crust and the lithosphere, we use the lithospheric thickness from {cite:t}`priestley:etal:2018`
and crustal thicknesses from the crust1.0 model {cite}`laske:elal:2012`.

Below "uppermost mantle thickness", the temperature and the density fields are computed using the corresponding depth-dependent scaling
values taken from {cite}`steinberger:calderwood:2006` relative to the input Vs anomalies (compositional field `Vs_anomaly`), using the
scaling files `dT_vs_scaling.txt` and `rho_vs_scaling.txt`, respectively, located in `cookbooks/tomography_based_plate_motions/input_data/`.

In order to avoid jumps in the temperature distribution, we smooth the temperatures between the two models across the "uppermost mantle thickness"
depth using a sigmoid function with a half-width of 20km.

The viscosity in the model follows a dislocation diffusion creep law with prefactors, activation energies and volumes for major phase
transitions. Additionally, we scale the laterally averaged viscosity to a reference viscosity profile from the profile of
{cite:t}`steinberger:calderwood:2006`. The user can choose from the available reference viscosity profiles in the
`cookbooks/tomography_based_plate_motions/input_data/viscosity_profiles` folder. The model defines different viscosities at the plate
boundaries and within the cratonic regions. For Earth-like behavior, plate boundaries and cratons have reduced and increased viscosity,
respectively, compared to the surrounding lithosphere.

## Running the model
The plugin that implements the material properties described above is `tomography_based_plate_motions.cc` located in the
[cookbooks/tomography_based_plate_motions/plugins](https://www.github.com/geodynamics/aspect/blob/main/cookbooks/tomography_based_plate_motions/plugins)
directory, and can be compiled with `cmake . && make` in this directory.
The model can then be run using the prm file [2D_slice_with_faults_and_cratons.prm](https://www.github.com/geodynamics/aspect/blob/main/cookbooks/tomography_based_plate_motions/2D_slice_with_faults_and_cratons.prm).
The sections relevant to this cookbook are:

```{literalinclude} tomography_based_plate_motions.part.prm
```

After running this model, we can visualize the computed material properties and analyze the instantaneous mantle flow.
{numref}`fig:mantle_flow` shows the computed non-adiabatic density and viscosity fields in the model along with the velocities generated
from the force balance between driving and resisting forces.

```{figure-md} fig:mantle_flow
<img src="Fig_mantle_flow.*" />

Instantaneous velocity field in the model overlying the heterogeneous lateral and radial viscosity distribution (left) and the non-adiabatic
densities (right). The magenta and white contours represent the compositional field corresponding to plate boundaries and cratons, respectively.
The cyan contour marks the lithospheric depths in our model, which also controls the depth until which the viscosities at the plate boundaries extend.
```
It is difficult to compare the modeled surface deformation with the observations using a 2D shell, but we can already see the expected
flow patterns such as upwellings due the positive buoyancy beneath the African plate and downwellings due to slabs under the Pacific plate.

The user can change several parameters to investigate the mantle forces and how they affect the flow, such as:

- The input tomography model, several options are available on the [IRIS DMC webpage](http://ds.iris.edu/ds/products/emc-earthmodels/).
- The density or temperature scaling relative to the input shear wave anomalies.
- The reference viscosity profile.
- The prescribed viscosity of the plate boundaries, which controls friction between plates.
- To use viscous and neutrally buoyant cratons in the model.

## Added complexity—Slab2 and initial topography
For additional complexity, we use the Slab2 {cite}`hayes2018slab2` database that describes in detail the three dimensional geometries
of all the seismically active subduction zones on Earth instead of the vertical slabs defined in the TM1 model. We also deform the outer
shell by the initial topography. Incorporating these finer heterogeneities makes the simulations more physically realistic and may improve
the fit to observed surface deformation patterns.

The slabs are included as a compositional field that is diffused with a length scale of 30 km. Physically, this value represents the
diffusion of a slab in $\approx$ 15 Ma using diffusivity of $1.14\times10^{-6} m^2 s^{-1}$. The diffusion not only approximates the thermal diffusion
of the slab as it warms up by the surrounding mantle, but also helps the solver convergence by avoiding sharp transition in material
properties between slabs and the surrounding mantle.
Within the slabs, we compute the temperatures using the diffused slabs' composition and a user-defined temperature anomaly, added to the
reference mantle adiabat. Additionally, we add a low-viscosity layer in the material model, the mid-mantle layer, that extends below the mantle
transition zone until 1000 km depth. The weak mid-mantle layer ensures that our slabs do not get stuck and can still pull the attached
lithospheric plates as they sink into the more viscous lower mantle.

For topography, we use a 1-deg spacing longitude-latitude grid and then smoothen the topography using a gaussian filter with standard deviation of 5.
The chosen smoothening reduces the pressure oscillations in our models that arise from the sharp topography gradients.

:::{note}
Our tests with initial topography revealed that having smooth topography variations is important for both linear and nonlinear
solver convergence. One possible explanation is that as the outer shell deforms according to the imposed topography, it is more
difficult to apply the free-slip boundary condition on a mesh that has sudden variations in the tangential vector across the
neighboring cells (see {numref}`fig:topography_test`).
:::

```{figure-md} fig:topography_test
<img src="Fig_topography_test.*" />

(a) Flow in a simple compressible spherical shell with added topography—a mountain and a basin. The arrows marks the region where the
tangential flow changes from deformed to undeformed surface. Shorter-wavelength topography changes leads to higher linear iterations
and higher non-linear residual in the system. (b) We smooth the observed topography in this cookbook because the material model is quite
complex already for sufficient solver convergence.
```

## Running the model
The prm file that includes all the heterogeneites is available in the folder,
[2D_slice_with_faults_slabs_and_topo.prm](https://www.github.com/geodynamics/aspect/blob/main/cookbooks/tomography_based_plate_motions/2D_slice_with_faults_slabs_and_topo.prm).
The relevant sections are :

```{literalinclude} tomography_based_plate_motions-2.part.prm
```

The parameter, `Use asthenosphere viscosity scaling in cold regions` ensures that the slabs are uniformly strong and that the low
viscosity of the asthenosphere does not affect the slabs (i.e, cold regions) during lateral viscosity averaging.

We can visualize the flow field around the Chilean slab in the {numref}`fig:flow_around_slabs` below:

```{figure-md} fig:flow_around_slabs
<img src="Fig_flow_around_slabs.*" />

Instantaneous velocity field illustrating subduction of the Chilean slab, slab pull on the attached Nazca plate, and resulting the
trench suction acting on the overriding South American plate. The low-viscosity layer added above the slab is an approximation of the
weak sediment layer above a slab. The white line represents the base of the low mid-mantle viscosity that allows the slabs to subduct
freely below the upper mantle.
```

## Surface deformation: plate motions
With access to high-performance computing, we would want to run these models in three dimensions and compare the model output to
surface observations such as plate motions or SHmax directions.
The following model ({numref}`fig:high_res_surface_deformation`) was run using ~2.5k cores on the NSF-funded High Performance
Computing System (Frontera) at the University of Texas, which shows a practical application of how these models can be used to approximate
the present-day physical state of the Earth.

```{figure-md} fig:high_res_surface_deformation
<img src="Fig_high_res_surface_deformation.*" />

Model setup (left) including all the heterogeneites and the modeled instantaneous plate velocities (right). The red and the black arrows
represent the modeled and the observed plate velocities, respectively.
```

> All the individual heterogeneous models are added in a modular fashion such that the user can easily modify whether or how to include them.
