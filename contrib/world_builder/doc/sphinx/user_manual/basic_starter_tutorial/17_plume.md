(part:user_manual:chap:basic_starter_tutorial:sec:17_plume)=
Adding a mantle plume
===============================

The last feature we will be adding in this tutorial is the `plume`. The plume is defined by a set of coordinates which represent where the plume is at different depths. The depths are set through the `cross section depths` parameter which requires as many depths as there are coordinates. Each cross-section is an ellipse for which the parameters can be set individually. Please see [the plume feature description](part:user_manual:chap:parameter_documentation:sec:features:subsec:plume) for more information on what each of the parameters do. In this case we will increase the eccentricity of the plume at the top.

With the definition of the area the feature plume contains completed, we can now add a temperature and compositional structure. A gaussian is a good first order approximation of the temperature, so we will add that as the temperature model. This temperature model allows you to set change the gaussian distribution parameters for each depth segment. Note that in this case we are not replacing the temperature, but adding to the temperature which was already there.


::::::{tab-set}

:::::{tab-item} Important lines
:sync: Partial

```{literalinclude} ../../_static/gwb_input_files/BST_17_plume.wb
:language: json
:lineno-start: 68
:lines: 68-88
:emphasize-lines: 4,5,6,7,11,12,13,14
```
::::{grid} 3
:::{grid-item-card} BST_17_plume.wb
:link: ../../_static/gwb_input_files/BST_17_plume.wb
:::
:::{grid-item-card} BST_17_plume.grid
:link: ../../_static/gwb_input_files/BST_17_plume.grid
:::
:::{grid-item-card} Paraview V4 state file 
:link: ../../_static/paraview_state_files/BST_v4.pvsm
:::
::::
:::::

:::::{tab-item} Full file
:sync: Full


```{literalinclude} ../../_static/gwb_input_files/BST_17_plume.wb
:language: json
:lineno-start: 1
:emphasize-lines: 71,72,73,74,78,79,80,81
```

::::{grid} 3
:::{grid-item-card} BST_17_plume.wb
:link: ../../_static/gwb_input_files/BST_17_plume.wb
:::
:::{grid-item-card} BST_17_plume.grid
:link: ../../_static/gwb_input_files/BST_17_plume.grid
:::
:::{grid-item-card} Paraview V4 state file 
:link: ../../_static/paraview_state_files/BST_v4.pvsm
:::
::::
:::::

::::::


```{figure} ../../../../doc/sphinx/_static/images/user_manual/basic_starter_tutorial/BST_17.png
:name: BST_17
:alt: Basic Starter Tutorial section 17 highres result. 
:align: center

Basic Starter Tutorial section 17 high resolution result, where the plume feature is used. In both the top and bottom figure, the area where the plume feature is present is shown. On the top this area is colored red, and on the bottom figure it is colored by temperature. This has 8 times the resolution then the grid file above.
```


This covers the full complexity of the World Builder model we are building for this tutorial. Well done for making it this far! You should now be able to start building your own models in the World Builder. In the next two tutorials, we are going to take a look at how to create a 2D model from this 3D model and how to make this into a spherical model.
