(part:user_manual:chap:basic_starter_tutorial:sec:07_subducting_plate_oceanic_part)=
Subducting plate oceanic part
========================

Now that we have made the overriding (Caribbean) plate, it is time to add the oceanic part of the subducting plate (the Atlantic). To do this, we just need to add another object to the feature list. For this plate we will assume the oceanic lithosphere is really old, and the temperature gradient in the lithosphere is therefore linear between a top temperature of 293.15K and the adiabatic temperature at the max depth of the plate. Luckily for us, this is the default in the World Builder, so we will only need to provide the `model` (`linear`) and the `max depth`, which we will set to 100km.

For the composition, we are going to do something a bit more special. A common thing you will probably want to do with compositional fields is have multiple layers of them within a feature. This is very easy to do in the GWB. If you remember from the last section, both the `temperature models` and `compositional models` are a list of objects. Making layers is thus as easy as adding multiple compositional models with each their own `min depth` and `max depth`. 

```{note}
If you provide a `min depth` or `max depth` outside the min and max depth range of the feature, it will be cut off by the feature min and max depth range.
```

::::::{tab-set}

:::::{tab-item} Important lines
:sync: Partial

```{literalinclude} ../../_static/gwb_input_files/BST_07_subducting_plate_oceanic_part.wb
:language: json
:lineno-start: 18
:lines: 18-25
:emphasize-lines: 5,6,7
```
::::{grid} 3
:::{grid-item-card} BST_07_subducting_plate_oceanic_part.wb
:link: ../../_static/gwb_input_files/BST_07_subducting_plate_oceanic_part.wb
:::
:::{grid-item-card} BST_07_subducting_plate_oceanic_part.grid
:link: ../../_static/gwb_input_files/BST_07_subducting_plate_oceanic_part.grid
:::
:::{grid-item-card} Paraview v2 state file 
:link: ../../_static/paraview_state_files/BST_v2.pvsm
:::
::::
:::::

:::::{tab-item} Full file
:sync: Full


```{literalinclude} ../../_static/gwb_input_files/BST_07_subducting_plate_oceanic_part.wb
:language: json
:lineno-start: 1
:emphasize-lines: 22,23,24
```

::::{grid} 3
:::{grid-item-card} BST_07_subducting_plate_oceanic_part.wb
:link: ../../_static/gwb_input_files/BST_07_subducting_plate_oceanic_part.wb
:::
:::{grid-item-card} BST_07_subducting_plate_oceanic_part.grid
:link: ../../_static/gwb_input_files/BST_07_subducting_plate_oceanic_part.grid
:::
:::{grid-item-card} Paraview v2 state file 
:link: ../../_static/paraview_state_files/BST_v2.pvsm
:::
::::
:::::

::::::


```{figure} ../../../../doc/sphinx/_static/images/user_manual/basic_starter_tutorial/BST_07.png
:name: BST_07
:alt: Basic Starter Tutorial section 7. 
:align: center

Basic Starter Tutorial section 7. The top part of the figure shows where the composition as been assigned as an object. Currently is shows composition 0 as green, composition 1 as yellow and composition 3 as blue. The bottom part shows the temperature as seen slightly from below where only temperatures between 300K and 1600K are shown. Both the oceanic plate with the ridge and the oceanic plate with a linear temperature profile are visible.
```
