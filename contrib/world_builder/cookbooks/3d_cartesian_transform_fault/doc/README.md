(part:user_manual:chap:cookbooks:sec:3d_cartesian_transform_fault)=
Mid-ocean ridge with a transform fault
=========================

The goal of this tutorial is to learn how to set up a model of a mid-ocean ridge with a transform fault, reproducing the setup of Behn et al., 2007: Thermal structure of oceanic transform faults. In addition, this tutorial shows how the created Geodynamic World Builder file can be used as initial condition for the geodynamic modeling code ASPECT. The corresponding cookbook for setting up the ASPECT model can be found [here](https://aspect-documentation.readthedocs.io/en/latest/user/cookbooks/cookbooks/transform_fault_behn_2007/doc/transform_fault_behn_2007.html).

The first step is to prescribe some global parameters for the GWB. This is important because we use the half-space cooling model to compute the temperature distribution within the oceanic plates, and this model uses the thermal diffusivity. In addition, the GWB needs to know what the adiabatic gradient is to compute the correct temperature, and we need to specify the surface and mantle temperature to be consistent with the study we want to reproduce.

Specifically, we here want to set the `surface temperature` to 0 degrees Celsius (273.25 K) and the `potential mantle temperature` to 1300 degrees Celsius (1573.25 K) as in Behn et al., 2007. Since that study does not include adiabatic heating, we need to set the `thermal expansion coefficient` to zero. To reproduce the temperature profile shown in Figure 2 in  Behn et al. (2007), we set the `thermal diffusivity` to 1.06060606e-6 m<sup>2</sup>/s (which corresponds to a thermal conductivity of 3.5 W/m/K assuming that the density is 3300 kg/m<sup>3</sup> and the specific keat is 1000 J/kg/K). Note that we need to make sure that these properties are consistent between the Geodynamic World Builder input file and the input file of the geodynamic modeling code that we are using this file as an initial condition for, in this case ASPECT.

::::::{tab-set}

:::::{tab-item} Global parameters
:sync: Partial

```{literalinclude} ../3d_cartesian_transform_fault.wb
:language: json
:lineno-start: 3
:lines: 3-6
```
::::{grid} 2
:::{grid-item-card} 3d_cartesian_transform_fault.wb
:link: ../3d_cartesian_transform_fault.wb
:::
:::{grid-item-card} 3d_cartesian_transform_fault.grid
:link: ../3d_cartesian_transform_fault.grid
:::
::::
:::::

:::::{tab-item} Full file
:sync: Full


```{literalinclude} ../3d_cartesian_transform_fault.wb
:language: json
:lineno-start: 1
:emphasize-lines: 3-6
```

::::{grid} 2
:::{grid-item-card} 3d_cartesian_transform_fault.wb
:link: ../3d_cartesian_transform_fault.wb
:::
:::{grid-item-card} 3d_cartesian_transform_fault.grid
:link: ../3d_cartesian_transform_fault.grid
:::
::::
:::::

::::::


As the next step, we prescribe the geometry and other properties of the oceanic plate. To do this, we create a feature of type `oceanic plate`, which we call `oceanic plate A`. We specify the `coordinates` where this plate should be present in the model. In this case, we want the whole model domain to be within the plate, so we choose coordinates that span the whole model domain, starting with the x=0 and y=0 corner and then going counterclockwise (increasing the x coordinate first) to all other corners at the surface. 

In addition to the location, we also need to specify the temperature distribution within the plate. 
To be able to do this, we need to provide the following information to the model: `model`, `ridge coordinates`, `spreading velocity` and `max depth`. Please see {ref}`open_features_items_oneOf_4_temperature-models_items_oneOf_2` for more info.
Here, we use a half-space cooling model (using the material properties we have already prescribed globally) as the temperature model. The full spreading rate is 6 cm/yr, so we need to prescribe a spreading velocity of 0.03 m/yr for the oceanic plate feature since the [spreading velocity](open_features_items_oneOf_4_temperature-models_items_oneOf_2_spreading-velocity) describes how fast the plate moves away from the ridge.
We also set the [top temperature](open_features_items_oneOf_4_temperature-models_items_oneOf_2_top-temperature) to 0 degrees Celsius (273.25 K), determining the temperature we use as the surface temperature in the half-space cooling model. Finally, we need to specify the coordinates of the mid-ocean ridge segments. Both segments are parallel to the y-axis, with the first segment being located at x=200 km, going from y=-1 km to y=50 km, and the second located at x=50 km going from y=50 km to y=101 km.


::::::{tab-set}

:::::{tab-item} Thermal structure
:sync: Partial

```{literalinclude} ../3d_cartesian_transform_fault.wb
:language: json
:lineno-start: 7
:lines: 7-24
:emphasize-lines: 4-14
```
::::{grid} 2
:::{grid-item-card} 3d_cartesian_transform_fault.wb
:link: ../3d_cartesian_transform_fault.wb
:::
:::{grid-item-card} 3d_cartesian_transform_fault.grid
:link: ../3d_cartesian_transform_fault.grid
:::
::::
:::::

:::::{tab-item} Full file
:sync: Full


```{literalinclude} ../3d_cartesian_transform_fault.wb
:language: json
:lineno-start: 1
:emphasize-lines: 10-20
```

::::{grid} 2
:::{grid-item-card} 3d_cartesian_transform_fault.wb
:link: ../3d_cartesian_transform_fault.wb
:::
:::{grid-item-card} 3d_cartesian_transform_fault.grid
:link: ../3d_cartesian_transform_fault.grid
:::
::::
:::::

::::::


The figures below show the resulting thermal structure of the model. 

```{figure} ./temperature_profile.png
:name: transform_profile
:alt: Transform fault temperature profile. 
:align: center

Temperature profile at the center of the transform fault. Output from the Geodynamic World Builder (orange dashed line) plotted on top of Figure 2a of Behn et al., 2007, showing good agreement with the half-space cooling model from that study.
```

```{figure} ./temperature_distribution.png
:name: transform_temperature
:alt: Transform fault temperature distribution. 
:align: center

Transform fault geometry and temperature distribution.
```
