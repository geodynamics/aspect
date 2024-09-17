(part:user_manual:chap:basic_starter_tutorial:sec:06_oceanic_plate_temperature)=
Oceanic plate temperature
=========================

Using uniform temperatures will only allow you to build very simple models. So let's make the model a bit more interesting by adding a half space cooling model to the overriding plate. This allows us to add a temperature structure based on the distance along the plate from a ridge. 

To be able to do this, we need to provide the following information to the model: `model`, `ridge coordinates`, `spreading velocity` and `max depth`. Please see {ref}`open_features_items_oneOf_4_temperature-models_items_oneOf_2` for more info.

We are already familiar with `max depth`, but `spreading velocity` and `ridge coordinates` are new. The [spreading velocity](open_features_items_oneOf_4_temperature-models_items_oneOf_2_spreading-velocity) is simply the velocity of how fast one part of the plate is moving away from the ridge. The [ridge coordinates](open_features_items_oneOf_4_temperature-models_items_oneOf_2_ridge-coordinates) are a list of spreading ridges, where each ridge is a list of 2D points. The half space is computed with the distance to the closest point on any ridge and the spreading velocity. 

If you are familiar with the half space cooling model, you would also think some kind of top and bottom/mantle temperature is required. These values can be provided, but by default they are set to reasonable values. The top temperature is set to 293.15K and the bottom temperature is set to -1K. Yes, that is negative 1 Kelvin! For this model, if the temperature is set to a negative value, it automatically uses the computed adiabatic temperature at that depth. This will be good enough for our tutorial, so we won't be changing those values here, but you are free to do so!


::::::{tab-set}

:::::{tab-item} Important lines
:sync: Partial

```{literalinclude} ../../_static/gwb_input_files/BST_06_oceanic_plate_temperature.wb
:language: json
:lineno-start: 6
:lines: 6-14
:emphasize-lines: 6,7
```
::::{grid} 3
:::{grid-item-card} BST_06_oceanic_plate_temperature.wb
:link: ../../_static/gwb_input_files/BST_06_oceanic_plate_temperature.wb
:::
:::{grid-item-card} BST_06_oceanic_plate_temperature.grid
:link: ../../_static/gwb_input_files/BST_06_oceanic_plate_temperature.grid
:::
:::{grid-item-card} Paraview v2 state file 
:link: ../../_static/paraview_state_files/BST_v2.pvsm
:::
::::
:::::

:::::{tab-item} Full file
:sync: Full


```{literalinclude} ../../_static/gwb_input_files/BST_06_oceanic_plate_temperature.wb
:language: json
:lineno-start: 1
:emphasize-lines: 11,12
```

::::{grid} 3
:::{grid-item-card} BST_06_oceanic_plate_temperature.wb
:link: ../../_static/gwb_input_files/BST_06_oceanic_plate_temperature.wb
:::
:::{grid-item-card} BST_06_oceanic_plate_temperature.grid
:link: ../../_static/gwb_input_files/BST_06_oceanic_plate_temperature.grid
:::
:::{grid-item-card} Paraview v2 state file 
:link: ../../_static/paraview_state_files/BST_v2.pvsm
:::
::::
:::::

::::::


```{figure} ../../../../doc/sphinx/_static/images/user_manual/basic_starter_tutorial/BST_06.png
:name: BST_06
:alt: Basic Starter Tutorial section 6. 
:align: center

Basic Starter Tutorial section 6. The top part of the figure shows where the composition has been assigned as an object. Currently it only shows composition 0 as green and is now limited to 100km depth. The bottom part shows the temperature as seen slightly from below where only temperatures between 300K and 1600K are shown. This allows a better view of the ridge.
```
