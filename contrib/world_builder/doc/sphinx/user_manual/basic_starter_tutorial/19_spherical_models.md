(part:user_manual:chap:basic_starter_tutorial:sec:19_spherical_models)=
Spherical models
================

We have already mentioned spherical coordinates in {ref}`part:user_manual:chap:basic_starter_tutorial:sec:03_coordinate_system`. In this section we will focus on converting our current Cartesian model to a spherical model. 

To convert, two changes need to be made:

1. The coordinate system should be set to spherical and a depth method needs to be set.
2. The coordinates you provide need to be changed from Cartesian to spherical.

To simplify the transition for this tutorial, we will just divide each coordinate by 10000 which will put this into a nice degree range.

::::::{tab-set}

:::::{tab-item} Important lines
:sync: Partial

```{literalinclude} ../../_static/gwb_input_files/BST_19_spherical_models.wb
:language: json
:lineno-start: 1
:lines: 1-10
:emphasize-lines: 3,9
```
::::{grid} 3
:::{grid-item-card} BST_19_spherical_models.wb
:link: ../../_static/gwb_input_files/BST_19_spherical_models.wb
:::
:::{grid-item-card} BST_19_spherical_models.grid
:link: ../../_static/gwb_input_files/BST_19_spherical_models.grid
:::
:::{grid-item-card} Paraview Spherical state file 
:link: ../../_static/paraview_state_files/BST_spherical.pvsm
:::
::::
:::::

:::::{tab-item} Full file
:sync: Full


```{literalinclude} ../../_static/gwb_input_files/BST_19_spherical_models.wb
:language: json
:lineno-start: 1
:emphasize-lines: 3,9,14,24,27,32,34,39,46
```

::::{grid} 3
:::{grid-item-card} BST_19_spherical_models.wb
:link: ../../_static/gwb_input_files/BST_19_spherical_models.wb
:::
:::{grid-item-card} BST_19_spherical_models.grid
:link: ../../_static/gwb_input_files/BST_19_spherical_models.grid
:::
:::{grid-item-card} Paraview Spherical state file 
:link: ../../_static/paraview_state_files/BST_spherical.pvsm
:::
::::
:::::

::::::


```{figure} ../../../../doc/sphinx/_static/images/user_manual/basic_starter_tutorial/BST_19.png
:name: BST_19
:alt: Basic Starter Tutorial section 19 highres result. 
:align: center

Basic Starter Tutorial section 19 high resolution result. This has 8 times the resolution than the grid file above. Note that you can see the Earth's curvature! 
```
