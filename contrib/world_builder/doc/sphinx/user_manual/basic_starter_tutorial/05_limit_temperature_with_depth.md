(part:user_manual:chap:basic_starter_tutorial:sec:05_limit_temperature_with_depth)=
Limit temperature with depth
============================


As you may have noticed, the temperature and compositional field 0 in the model both extend to the bottom of the model. However, our model is 500km deep and the oceanic lithosphere should only extend to 100km while the "crustal composition" should only extend to 50km in depth. Luckily, this is easily solved. We start by setting a global maximum depth of 100km for the whole feature by adding `"max depth":100e3` to the end of line 7. We can also set maximums for individual objects so let's limit the uniform composition by adding `"max depth":50e3` to line 12. Please note that units are in meters.

```{note}
If the `max depth` of the composition is larger than the feature `max depth`, the composition will be be cut off as it is limited by the feature `max depth`.
```

::::::{tab-set}

:::::{tab-item} Important lines
:sync: Partial

```{literalinclude} ../../_static/gwb_input_files/BST_05_limit_depth.wb
:language: json
:lineno-start: 5
:lines: 5-15
:emphasize-lines: 3,8
```
::::{grid} 3
:::{grid-item-card} BST_05_limit_depth.wb
:link: ../../_static/gwb_input_files/BST_05_limit_depth.wb
:::
:::{grid-item-card} BST_05_limit_depth.grid
:link: ../../_static/gwb_input_files/BST_05_limit_depth.grid
:::
:::{grid-item-card} Paraview v1 state file 
:link: ../../_static/paraview_state_files/BST_v1.pvsm
:::
::::
:::::

:::::{tab-item} Full file
:sync: Full


```{literalinclude} ../../_static/gwb_input_files/BST_05_limit_depth.wb
:language: json
:lineno-start: 1
:emphasize-lines: 7,12
```

::::{grid} 3
:::{grid-item-card} BST_05_limit_depth.wb
:link: ../../_static/gwb_input_files/BST_05_limit_depth.wb
:::
:::{grid-item-card} BST_05_limit_depth.grid
:link: ../../_static/gwb_input_files/BST_05_limit_depth.grid
:::
:::{grid-item-card} Paraview v1 state file 
:link: ../../_static/paraview_state_files/BST_v1.pvsm
:::
::::
:::::

::::::


```{figure} ../../../../doc/sphinx/_static/images/user_manual/basic_starter_tutorial/BST_05.png
:name: BST_05
:alt: Basic Starter Tutorial section 5. 
:align: center

Basic Starter Tutorial section 5. The top part of the figure shows where the composition has been assigned as an object. Composition 0 (green) is now limited to 50km in depth. The bottom part shows the temperature as viewed slightly from below. The temperature of the overriding plate (blue), 293K, is now limited to 100km in depth.
```
