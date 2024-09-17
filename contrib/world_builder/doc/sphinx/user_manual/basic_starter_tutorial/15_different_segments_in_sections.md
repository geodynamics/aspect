(part:user_manual:chap:basic_starter_tutorial:sec:15_different_segment_in_sections)=
Different segments in sections
===============================

Slabs and faults may vary in depth as well as laterally. To accomplish this in the World Builder, we will need to introduce the term "sections". Sections are a way to overwrite the segments at a specified coordinate. 

```{note}
One way of thinking about segments and sections is that together they form a 2D grid. The segments form one axis (down dip) and the sections form the other axis (along strike). 
```

**Sections replace all the segments at that coordinate. It is required that all sections have the same number of segments**.

In this case we are going to change the values for the first coordinate, which is coordinate 0. We keep the length the same as the default first segment (compare lines 48 and 75), but we make the second segment a bit shorter (compare lines 52 and 58). We also do not change the feature composition model in the first segment, but we do change it in the second segment where we set it equal to composition of the lower part of the oceanic plate for the whole thickness of the segment (line 59).

::::::{tab-set}

:::::{tab-item} Important lines
:sync: Partial

```{literalinclude} ../../_static/gwb_input_files/BST_15_different_segments_in_sections.wb
:language: json
:lineno-start: 43
:lines: 43-63
:emphasize-lines: 12,13,14,15,16,17,18
```
::::{grid} 3
:::{grid-item-card} BST_15_different_segments_in_sections.wb
:link: ../../_static/gwb_input_files/BST_15_different_segments_in_sections.wb
:::
:::{grid-item-card} BST_15_different_segments_in_sections.grid
:link: ../../_static/gwb_input_files/BST_15_different_segments_in_sections.grid
:::
:::{grid-item-card} Paraview v3 state file 
:link: ../../_static/paraview_state_files/BST_v3.pvsm
:::
::::
:::::

:::::{tab-item} Full file
:sync: Full


```{literalinclude} ../../_static/gwb_input_files/BST_15_different_segments_in_sections.wb
:language: json
:lineno-start: 1
:emphasize-lines: 54,55,56,57,58,59,60
```

::::{grid} 3
:::{grid-item-card} BST_15_different_segments_in_sections.wb
:link: ../../_static/gwb_input_files/BST_15_different_segments_in_sections.wb
:::
:::{grid-item-card} BST_15_different_segments_in_sections.grid
:link: ../../_static/gwb_input_files/BST_15_different_segments_in_sections.grid
:::
:::{grid-item-card} Paraview v3 state file 
:link: ../../_static/paraview_state_files/BST_v3.pvsm
:::
::::
:::::

::::::

```{figure} ../../../../doc/sphinx/_static/images/user_manual/basic_starter_tutorial/BST_15.png
:name: BST_15
:alt: Basic Starter Tutorial section 15 highres result. 
:align: center

Basic Starter Tutorial section 15 high resolution result. This has 8 times the resolution then the grid file above. Note that the part of the slab on the far side now has a different composition and angle.
```
