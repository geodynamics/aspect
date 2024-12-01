(part:user_manual:chap:basic_starter_tutorial:sec:13_subducting_slab_adding_a_segment)=
Subducting slab: adding a segment
================================

Although our slabs can now bend, real slabs actually can bend many times in many different ways. The way this is parameterized for faults and slab in the World Builder is through dividing the slab into multiple segments. A segment is a downwards section of a slab or a fault. Each segment can have their own parameters and models attached to it. 

So how do we achieve this? We just add a second object to the segments line. Usually you will want to make sure that the end of one segment has the same angle and thickness as the beginning of the next segment. So in this case we end the first segment with an angle of 60 degrees and begin the second segment also with an angle of 60 degrees.

::::::{tab-set}

:::::{tab-item} Important lines
:sync: Partial

```{literalinclude} ../../_static/gwb_input_files/BST_13_subducting_slab_adding_a_segment.wb
:language: json
:lineno-start: 43
:lines: 43-52
:emphasize-lines: 7
```
::::{grid} 3
:::{grid-item-card} BST_13_subducting_slab_adding_a_segment.wb
:link: ../../_static/gwb_input_files/BST_13_subducting_slab_adding_a_segment.wb
:::
:::{grid-item-card} BST_13_subducting_slab_adding_a_segment.grid
:link: ../../_static/gwb_input_files/BST_13_subducting_slab_adding_a_segment.grid
:::
:::{grid-item-card} Paraview v3 state file 
:link: ../../_static/paraview_state_files/BST_v3.pvsm
:::
::::
:::::

:::::{tab-item} Full file
:sync: Full


```{literalinclude} ../../_static/gwb_input_files/BST_13_subducting_slab_adding_a_segment.wb
:language: json
:lineno-start: 1
:emphasize-lines: 49
```

::::{grid} 3
:::{grid-item-card} BST_13_subducting_slab_adding_a_segment.wb
:link: ../../_static/gwb_input_files/BST_13_subducting_slab_adding_a_segment.wb
:::
:::{grid-item-card} BST_13_subducting_slab_adding_a_segment.grid
:link: ../../_static/gwb_input_files/BST_13_subducting_slab_adding_a_segment.grid
:::
:::{grid-item-card} Paraview v3 state file 
:link: ../../_static/paraview_state_files/BST_v3.pvsm
:::
::::
:::::

::::::

```{figure} ../../../../doc/sphinx/_static/images/user_manual/basic_starter_tutorial/BST_13.png
:name: BST_13
:alt: Basic Starter Tutorial section 13. 
:align: center

Basic Starter Tutorial section 13. The top part of the figure shows where the composition has been assigned as an object. Currently it shows composition 0 as green, composition 1 as yellow, composition 2 as purple and composition 3 as blue. Composition 4 is not shown to be able to see the slab. The front half of the overriding plate (composition 1) has also been removed to be able to better view the slab. The bottom part shows the temperature as seen slightly from below where only temperatures between 300K and 1600K are shown. 
```
