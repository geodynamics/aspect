(part:user_manual:chap:basic_starter_tutorial:sec:08_passive_margin_variable_depth)=
Adding a passive margin with variable depth
============================

Although we have hinted towards it before, there is also an [area feature](part:user_manual:chap:concepts:sec:area_features) called `continental plate`. It works as you would expect based on what you have seen with the oceanic plate. Here we are going to use it to do something very cool in our model - specify a feature with variable depths for the whole feature and also layers within the given feature!

To showcase this, we will be adding a passive margin to our model at the overriding plate side of our model. To start out, we need to know that both the `min depth` and `max depth` accept two types of values: a number and an array of values at points. 

```{note}
For an example of what this looks like in technical terms, see {ref}`open_features_items_oneOf_1_max-depth`.
```

The number input we have seen before and it sets the maximum depth of the feature to a single value. The value at points system works a bit differently and can be used to achieve the same result. To start, this is exactly what we are going to do for the `max depth` of the continental feature.

The value at points system creates a list of points that includes all the initial edge points of an object and assigns a value to them. We then interpolate between close points to obtain the actual values at the remaining points according to the model type, e.g., linearly to the points set in `max depth` below.

Below are two examples of how to first set the depth of two points to 200 km, and then an example of how to also set a third point to a value of 100 km.

```{code-block} json
---
lineno-start: 1
---
    "max depth":[[200e3, [[10e3,0],[20e3,10e3]]]],
    
    "max depth":[[200e3, [[10e3,0],[20e3,10e3]]],[100e3,[[15e3,0]]]]
```

A common operation is to set all the corners to a single value, and then maybe overwrite individual corners. To make life easy, there is a quick way to (re)set the value at the corner points: just pass a value with no points. This is what is done in the code sample below at the emphasized lines. The first entry is a value without points. That means that all the corner points are set to that value.

```{note}
If you provide a point twice, explicitly, or implicitly through the use of the corner values shortcut (e.g. `[200]`), the last defined value is used. This follows the painting analogy used before, where you overpaint older values.
```

::::::{tab-set}

:::::{tab-item} Important lines
:sync: Partial

```{literalinclude} ../../_static/gwb_input_files/BST_08_passive_margin_variable_depth.wb
:language: json
:lineno-start: 16
:lines: 16-30
:emphasize-lines: 2,6,11,13
```
::::{grid} 3
:::{grid-item-card} BST_08_passive_margin_variable_depth.wb
:link: ../../_static/gwb_input_files/BST_08_passive_margin_variable_depth.wb
:::
:::{grid-item-card} BST_08_passive_margin_variable_depth.grid
:link: ../../_static/gwb_input_files/BST_08_passive_margin_variable_depth.grid
:::
:::{grid-item-card} Paraview v2 state file 
:link: ../../_static/paraview_state_files/BST_v2.pvsm
:::
::::
:::::

:::::{tab-item} Full file
:sync: Full


```{literalinclude} ../../_static/gwb_input_files/BST_08_passive_margin_variable_depth.wb
:language: json
:lineno-start: 1
:emphasize-lines: 17,21,26,28
```

::::{grid} 3
:::{grid-item-card} BST_08_passive_margin_variable_depth.wb
:link: ../../_static/gwb_input_files/BST_08_passive_margin_variable_depth.wb
:::
:::{grid-item-card} BST_08_passive_margin_variable_depth.grid
:link: ../../_static/gwb_input_files/BST_08_passive_margin_variable_depth.grid
:::
:::{grid-item-card} Paraview v2 state file 
:link: ../../_static/paraview_state_files/BST_v2.pvsm
:::
::::
:::::

::::::


```{figure} ../../../../doc/sphinx/_static/images/user_manual/basic_starter_tutorial/BST_08.png
:name: BST_08
:alt: Basic Starter Tutorial section 8. 
:align: center

Basic Starter Tutorial section 8. The top part of the figure shows where the composition has been assigned as an object. Currently it shows composition 0 as green, composition 1 as yellow and composition 3 as blue. The bottom part shows the temperature as seen slightly from below where only temperatures between 300K and 1600K are shown. The added continental plate with variable thickness of its two layers is now visible on the left side of the image.
```
