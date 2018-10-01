# Mantle convection using tomography data

Mantle flow models based on present-day observational constraints can help us to investigate several global and regional tectonic processes.
This cookbook describes how to setup such a model based on input tomography data and plate boundary locations. These models are typically computationally very expensive in 3D to make meaningful comparisons with the observations, therefore, in this cookbook we explain the workflow for a 2D spherical shell.

## Description of input data

### Tomography file
In this cookbook, the input tomography data is LLNL-G3D-JPS{cite}`simmons:etal:2019`, available at the [IRIS DMC webpage](http://ds.iris.edu/ds/products/emc-earthmodels/), read in as an ASCII compositional data. The input file is located in `cookbooks/equilibirum_grain_size/input_data/` and has the following format

```{literalinclude} input_tomography.part.txt
```

where grain_size is set to a constant value of 5mm, Vp and Vs are the P- and S-wave seismic velocities, respectively, Vs_anomaly is the S-wave velocity anomaly relative a reference velocity profile, faults represent the composition value at plate boundaries, and cratons represent the composition value within the continental cratons. Both faults and cratons have composition values set to 0 in this ASCII file because they are defined separately using Worldbuilder{cite}`fraters_menno_2020_3900603` and the input file `cookbooks/equilibrium_grain_size/input_file/world_builder_smac_cratons_faults_2D.json`.

### World Builder file
We use the plate boundary geometry defined in the Nuvel model{cite}`demets:etal:1990` input as `fault` features and the locations of cratons are taken from{cite}`nataf:etal:1996`, input as `continental plate` features regions in WorldBuilder (see: reference for more details). The compositional value of faults transitions smoothly from 1 at the fault trace to 0 over a width of 50 km on either side of the fault. We do this to ensure that the prescribed plate boundary viscosity which is based on the fault composition is also smooth.
Note that while the input file includes plate boundary geometry and cratons for the whole Earth, we only compute the compositions along a cross-section in this cookbook. Specifically, we the coordinates in `cross-section` represents the compositions along the equator. The input faults and cratons in a 3D spherical model would look like {numref}`fig:composition`.

```{figure-md} fig:composition
<img src="Fig_composition.*" />

The composition of plate boundary geometry ("faults") and cratons ("continents") using the input file coordinates for a 3D model.
```

## Computed material properties
A common practice in mantle convection models using tomography data is to use a different model for the uppermost mantle structure because of the limited seimic resolution to accurately resolve the lithospheric structure and slabs. For this cookbook, we use a temperature model, TM1{cite}`osei:etal:2018`, located in `cookbooks/equilibrium_grain_size/input_data/upper_mantle_TM1_2D.txt`. The user can control the depth until which we want to use this temperature model using the `uppermost mantle thickness` parameter.
{numref}`fig:Vs_temperatures` shows the input Vs anomalies and the computed temperatures in the model.

```{figure-md} fig:Vs_temperatures
<img src="Fig_tomography_temperature.*" />

The input Vs anomalies (left) and the computed temperatures (right). Note that we only use tomography-based temperatures below the `uppermost mantle thickness`  depth, marked by the white contour.
```

We compute the density from the thermal anomalies in the TM1 model relative to a reference temperature of 293 K as done by {cite:t}`osei:etal:2018`. We use constant thermal expansion coefficients and comprehensibilities within the crust, lithospheric mantle and asthenosphere, respectively. To define the crust and the lithosphere, we use the lithospheric thickness from{cite}`priestley:etal:2018` and crustal thicknesses from the crust1.0 model{cite}`laske:elal:2012`.

Below `uppermost mantle thickness`, the temperature and the density fields are computed using the corresponding depth-dependent scaling values taken from{cite}`steinberger:calderwood:2006` relative to the input Vs anomalies (compositional field `Vs_anomaly`),using the scaling files `dT_vs_scaling.txt` and `rho_vs_scaling.txt`, respectively, located in `cookbooks/equilibirum_grain_size/input_data/`.

In order to avoid jumps in the temperature distribution, we smooth the temperatures between the two models across the `uppermost mantle thickness` depth.

The viscosity in the model follows a dislocation diffusion creep with prefactors, activation energies and volumes for major phase transitions. Additionally, we scale the laterally averaged viscosity to a reference viscosity profile, in particular, we use{cite:t}`steinberger:calderwood:2006`'s profile. The user can choose from the available reference viscosity profiles in the `cookbooks/equilibirum_grain_size/input_data/viscosity_profiles` folder.
The model defines different viscosities at the plate boundaries and  within the cratonic regions. For Earth-like behavior, plate boundaries and cratons have reduced and increased viscosity, respectively, compared to the surrounding lithosphere.

## Running the model
The plugin that implements the material properties described above can be found in [cookbooks/equilibrium_grain_size/plugins/equilibrium_grain_size.cc][] and can be compiled with
`cmake . && make` in the [cookbooks/equilibrium_grain_size/plugins/][] directory.
The model can then be run using the prm file [cookbooks/equilibrium_grain_size/models/tests/2D_slice_with_faults_and_cratons.prm](https://www.github.com/geodynamics/aspect/blob/main/cookbooks/equilibirum_grain_size/models/tests/2D_slice_with_faults_and_cratons.prm). The sections relevant to this cookbook are:

```{literalinclude} tomography_based_model.part.prm
```

After running this model, we can visualize the computed material properties and analyze the instantaneous mantle flow. {numref}`fig:mantle_flow` shows the computed non-adiabatic density and viscosity fields in the model along with the velocities generated from the force balance between driving and resisting forces.

```{figure-md} fig:mantle_flow
<img src="Fig_mantle_flow.*" />

Instantaneous velocity field in the model overlying the heterogeneous lateral and radial viscosity distribution (left) and the non-adiabatic densities (right). The magenta and white contours represent the compositional field corresponding to plate boundaries and cratons, respectively. Cyan contour mark the lithospheric depths in our model, which also controls the depth until which the viscosities at the plate boundaries extend.
```
It is difficult to make comparison of the surface deformation with the observations using a 2D shell, but we can already see the expected flow patterns such as upwellings due the positive buoyancy beneath the African plate and downwellings due to slabs under the Pacific plate.

The user can change several parameters to investigate the mantle forces and how they affect the flow, such as:

- The input tomography model, several options are available on the [IRIS DMC webpage](http://ds.iris.edu/ds/products/emc-earthmodels/).
- The density or temperature scaling relative to the input shear wave anomalies.
- The reference viscosity profile.
- The prescribed viscosity of the plate boundaries, which controls friction between plates.
- To use viscous and neutrally buoyant cratons in the model.
