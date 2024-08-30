(part:user_manual:chap:basic_starter_tutorial:sec:18_2D_models)=
2D models
=========

2D models in the World Builder are nothing more than a cross section through a 3D model. To create a 2D model, you must state in the World Builder file the origin of the cross section (where x=0 and depth=0) and in what direction the cross section should go (the positive x direction). You can set this with the global parameter `cross section`. It takes two points. The first point is the origin, and the second point is the direction of the cross section. 

```{note}
The origin of the 3D model doesn't have to be the same as the origin of the 2D cross section. In the model below, the origin of the 2D cross section is located at [0,450e3] (plane view, depth=0).
```

```{note}
Even if you have a `cross section` defined in your World Builder file, you can still use it for 3D models. If you want to use 2D models, you will need to have a cross section defined.
```

Now, let's add a cross section through the slab in our model. 

::::::{tab-set}

:::::{tab-item} Important lines
:sync: Partial

```{literalinclude} ../../_static/gwb_input_files/BST_18_2D_models.wb
:language: json
:lineno-start: 1
:lines: 1-6
:emphasize-lines: 4
```
::::{grid} 3
:::{grid-item-card} BST_18_2D_models.wb
:link: ../../_static/gwb_input_files/BST_18_2D_models.wb
:::
:::{grid-item-card} BST_18_2D_models.grid
:link: ../../_static/gwb_input_files/BST_18_2D_models.grid
:::
:::{grid-item-card} Paraview 2D state file 
:link: ../../_static/paraview_state_files/BST_2D.pvsm
:::
::::
:::::

:::::{tab-item} Full file
:sync: Full


```{literalinclude} ../../_static/gwb_input_files/BST_18_2D_models.wb
:language: json
:lineno-start: 1
:emphasize-lines: 4
```

::::{grid} 3
:::{grid-item-card} BST_18_2D_models.wb
:link: ../../_static/gwb_input_files/BST_18_2D_models.wb
:::
:::{grid-item-card} BST_18_2D_models.grid
:link: ../../_static/gwb_input_files/BST_18_2D_models.grid
:::
:::{grid-item-card} Paraview 2D state file 
:link: ../../_static/paraview_state_files/BST_2D.pvsm
:::
::::
:::::

::::::

```{note}
You need to change the gridfile to a 2D grid to be able to see the difference.
```


```{figure} ../../../../doc/sphinx/_static/images/user_manual/basic_starter_tutorial/BST_18_3D_cross_section.png
:name: BST_18_3D_cross_section
:alt: Basic Starter Tutorial section 18 highres result. 
:align: center

The location of the 2D cross-section in the 3D model. The turquoise arrow shows the cross section origin, the dot at `[0,450e3]`, and direction, arrow head at `[100e3,450e3]` of the 2D plane. The semi-transparent turquoise plane show the full cross section location.
```

```{figure} ../../../../doc/sphinx/_static/images/user_manual/basic_starter_tutorial/BST_18.png
:name: BST_18_highres
:alt: Basic Starter Tutorial section 18 highres result. 
:align: center

Basic Starter Tutorial section 18 high resolution result. This has 4 times the resolution than the grid file above. Note that some of the issues with the slab, like its abrubt ending can be solved by using the mass conserving temperature model instead of the McKenzie plate model. 
```
