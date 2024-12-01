(part:user_manual:chap:basic_starter_tutorial:sec:11_dip_change_in_segment)=
Dip change in segment
=====================

What we have achieved in the last section was already great but slabs in the Earth are usually not thought of as straight lines.  For example, near the surface they start with a dip angle of zero and then may increase their dip with depth to maybe 60 degrees. So, how do we do this in the World Builder?

This is where the array for the `angle` comes in. Last time we provided just one value but we are also allowed to provide two values. If you provide two values, the first value is the dip angle at the top of the segment and the second value is the dip angle at the bottom of the segment. To achieve this we only need to add a zero before the 60 in the `angle` array.


::::::{tab-set}

:::::{tab-item} Important lines
:sync: Partial

```{literalinclude} ../../_static/gwb_input_files/BST_11_dip_change_in_segment.wb
:language: json
:lineno-start: 43
:lines: 43-49
:emphasize-lines: 4
```
::::{grid} 3
:::{grid-item-card} BST_11_dip_change_in_segment.wb
:link: ../../_static/gwb_input_files/BST_11_dip_change_in_segment.wb
:::
:::{grid-item-card} BST_11_dip_change_in_segment.grid
:link: ../../_static/gwb_input_files/BST_11_dip_change_in_segment.grid
:::
:::{grid-item-card} Paraview v3 state file 
:link: ../../_static/paraview_state_files/BST_v3.pvsm
:::
::::
:::::

:::::{tab-item} Full file
:sync: Full


```{literalinclude} ../../_static/gwb_input_files/BST_11_dip_change_in_segment.wb
:language: json
:lineno-start: 1
:emphasize-lines: 46
```

::::{grid} 3
:::{grid-item-card} BST_11_dip_change_in_segment.wb
:link: ../../_static/gwb_input_files/BST_11_dip_change_in_segment.wb
:::
:::{grid-item-card} BST_11_dip_change_in_segment.grid
:link: ../../_static/gwb_input_files/BST_11_dip_change_in_segment.grid
:::
:::{grid-item-card} Paraview v3 state file 
:link: ../../_static/paraview_state_files/BST_v3.pvsm
:::
::::
:::::

::::::

```{note}
You can do the same for the `thickness`. If we want to start with a thickness of 100km and linearly taper it down to 0 km to form a slab tip, we can do the following: `"thickness":[100e3,0]`
```

```{todo}
The explanation can be significantly improved by adding conceptual figures
```

```{figure} ../../../../doc/sphinx/_static/images/user_manual/basic_starter_tutorial/BST_11.png
:name: BST_11
:alt: Basic Starter Tutorial section 11. 
:align: center

Basic Starter Tutorial section 11. The top part of the figure shows where the composition has been assigned as an object. Currently is shows composition 0 as green, composition 1 as yellow, composition 2 as purple and composition 3 as blue. Composition 4 is not shown to be able to see the slab. The front half of the overriding plate (composition 1) has also been removed to be able to better view the slab. The bottom part shows the temperature as seen slightly from below where only temperatures between 300K and 1600K are shown. 
```
