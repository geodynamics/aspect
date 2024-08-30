(part:user_manual:chap:basic_starter_tutorial:sec:14_different_models_in_segment)=
Different models in segments
============================

Slabs and faults can be heterogeneous in depth. This is very easy to represent in the World Builder. The segments we just added allow for defining temperature and compositional models that overwrite the temperature and compositional models for the whole feature. You can think of it this way (remember the painting analogy!): you set the feature models for each segment in the feature (line 54 and 55) and then you overwrite that feature model in the segments where you want something else. 

In this case we want to have two layers in the first segment where the upper layer is the same as the upper layer of the oceanic plate and the lower layer is the same as the rest of the subducting plate. This layout is mainly chosen to highlight how the transition between the oceanic plate can both be seamless and sharp.

::::::{tab-set}

:::::{tab-item} Important lines
:sync: Partial

```{literalinclude} ../../_static/gwb_input_files/BST_14_different_models_in_segments.wb
:language: json
:lineno-start: 43
:lines: 43-56
:emphasize-lines: 7,8,9
```
::::{grid} 3
:::{grid-item-card} BST_14_different_models_in_segments.wb
:link: ../../_static/gwb_input_files/BST_14_different_models_in_segments.wb
:::
:::{grid-item-card} BST_14_different_models_in_segments.grid
:link: ../../_static/gwb_input_files/BST_14_different_models_in_segments.grid
:::
:::{grid-item-card} Paraview v3 state file 
:link: ../../_static/paraview_state_files/BST_v3.pvsm
:::
::::
:::::

:::::{tab-item} Full file
:sync: Full


```{literalinclude} ../../_static/gwb_input_files/BST_14_different_models_in_segments.wb
:language: json
:lineno-start: 1
:emphasize-lines: 49,50,51
```

::::{grid} 3
:::{grid-item-card} BST_14_different_models_in_segments.wb
:link: ../../_static/gwb_input_files/BST_14_different_models_in_segments.wb
:::
:::{grid-item-card} BST_14_different_models_in_segments.grid
:link: ../../_static/gwb_input_files/BST_14_different_models_in_segments.grid
:::
:::{grid-item-card} Paraview v3 state file 
:link: ../../_static/paraview_state_files/BST_v3.pvsm
:::
::::
:::::

::::::


```{figure} ../../../../doc/sphinx/_static/images/user_manual/basic_starter_tutorial/BST_14.png
:name: BST_14
:alt: Basic Starter Tutorial section 14. 
:align: center

Basic Starter Tutorial section 14. The top part of the figure shows where the composition has been assigned as an object. Currently it shows composition 0 as green, composition 1 as yellow, composition 2 as purple and composition 3 as blue. Composition 4 is not shown to be able to see the slab. The front half of the overriding plate (composition 1) has also been removed to be able to better view the slab. The bottom part shows the temperature as seen slightly from below where only temperatures between 300K and 1600K are shown. 
```
