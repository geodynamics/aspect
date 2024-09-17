(part:user_manual:chap:basic_starter_tutorial:sec:09_optional_adding_mantle_layer)=
Optional adding mantle layer
============================

Sometimes you want precise control over what is happening in the mantle. In this case, we would want to have a specific composition for the mantle. 

```{note}
This is not always needed. For example in the visco-plastic material model in the FEM ASPECT, if no composition is defined, it is assumed to be a background composition. If you don't define a composition in an area in the World Builder, the World Builder will tell anyone who wants to know the composition at that location that the fraction of every composition is zero.
```

We are going to add the mantle layer as the first object in the feature list so that subsequent features can override it if desired. The model name for a mantle layer is `mantle layer`, and we will name it "upper mantle". Because we already have used compositions 0, 1 and 3 and will use 2 later, we will set this object to composition 4. We could also specify a different temperature model. The default is `adiabatic` but you can easily set it to something else or change the parameters of the adiabatic gradient.

```{note}
The [potential mantle temperature](open_features_items_oneOf_3_temperature-models_items_oneOf_1_potential-mantle-temperature), the [thermal expansion coefficient](open_features_items_oneOf_3_temperature-models_items_oneOf_1_thermal-expansion-coefficient) and the [specific heat](open_features_items_oneOf_3_temperature-models_items_oneOf_1_specific-heat) here have a default value of -1. This is, like we have seen before, a special value. In this case it means that the global value is used (as you can read in the documentation of the value). It is **strongly recommended to only set the global values** so that the whole setup remains self-consistent, but if you do want to change it in individual locations, you can!
```

::::::{tab-set}

:::::{tab-item} Important lines
:sync: Partial

```{literalinclude} ../../_static/gwb_input_files/BST_09_adding_mantle_layer.wb
:language: json
:lineno-start: 5
:lines: 5-10
:emphasize-lines: 5
```
::::{grid} 3
:::{grid-item-card} BST_09_adding_mantle_layer.wb
:link: ../../_static/gwb_input_files/BST_09_adding_mantle_layer.wb
:::
:::{grid-item-card} BST_09_adding_mantle_layer.grid
:link: ../../_static/gwb_input_files/BST_09_adding_mantle_layer.grid
:::
:::{grid-item-card} Paraview v2 state file 
:link: ../../_static/paraview_state_files/BST_v2.pvsm
:::
::::
:::::

:::::{tab-item} Full file
:sync: Full


```{literalinclude} ../../_static/gwb_input_files/BST_09_adding_mantle_layer.wb
:language: json
:lineno-start: 1
:emphasize-lines: 9
```

::::{grid} 3
:::{grid-item-card} BST_09_adding_mantle_layer.wb
:link: ../../_static/gwb_input_files/BST_09_adding_mantle_layer.wb
:::
:::{grid-item-card} BST_09_adding_mantle_layer.grid
:link: ../../_static/gwb_input_files/BST_09_adding_mantle_layer.grid
:::
:::{grid-item-card} Paraview v2 state file 
:link: ../../_static/paraview_state_files/BST_v2.pvsm
:::
::::
:::::

::::::

```{figure} ../../../../doc/sphinx/_static/images/user_manual/basic_starter_tutorial/BST_09.png
:name: BST_09
:alt: Basic Starter Tutorial section 9. 
:align: center

Basic Starter Tutorial section 9. The top part of the figure shows where the composition has been assigned as an object. Currently it shows composition 0 as green, composition 1 as yellow, composition 3 as blue and composition 4 is shown as red. The bottom part shows the temperature as seen slightly from below where only temperatures between 300K and 1600K are shown. We will not show composition 4 in future figures anymore to be able to show other features more clearly.
```
