(part:user_manual:chap:basic_starter_tutorial:sec:16_mass_conserving)=
Changing to a mass conserving slab temperature
===============================

In the section on [slab temperature](part:user_manual:chap:basic_starter_tutorial:sec:12_subducting_plate_temp), we added the {cite:t}`McKenzie_1970` (that is the plate model) slab temperature structure. Although this is a good first order approximation of a slab temperature, the recently develop {cite:t}`Billen_Fraters_AGU_2023` temperature model has many advantages over it, with the downside that is it a bit more involved to use. In this chapter we will just switch out the plate model with mass conserving without going too much into the detail of how to actually use it in practice. If you plan to use the mass conserving temperature model, please first read through both the [Simple Subduction Model: 2D Cartesian](part:user_manual:chap:cookbooks:sec:simple_subduction_2d_cartesian)and [Simple Subduction Model: 2D Chunk](part:user_manual:chap:cookbooks:sec:simple_subduction_2d_chunk) cookbooks.

## Changing the Subducting oceanic plate temperature
One of the advantages of using the mass conserving slab temperature model is that it doesn't assume a linear temperature structure at the trench. This means we can seamlessly connect a half space cooling model or a plate model to the mass conserving slab. But for that to make sense, we need to change the subucting oceanic plate to a half space cooling model (the default for the mass conserving temperature model). In this case we will put the ridge far away. Because the half space model can affect the temperature at much deeper depths, we also need to change the max depth for the feature and models. Note that we will want to keep the max depth of the composition at 100km, so we now need to set a max depth in the lowest layer of the composition model.


::::::{tab-set}

:::::{tab-item} Important lines
:sync: Partial

```{literalinclude} ../../_static/gwb_input_files/BST_16_mass_conserving.wb
:language: json
:lineno-start: 36
:lines: 36-43
:emphasize-lines: 2,4,5,7
```
::::{grid} 3
:::{grid-item-card} BST_16_mass_conserving.wb
:link: ../../_static/gwb_input_files/BST_16_mass_conserving.wb
:::
:::{grid-item-card} BST_16_mass_conserving.grid
:link: ../../_static/gwb_input_files/BST_16_mass_conserving.grid
:::
:::{grid-item-card} Paraview 2D state file 
:link: ../../_static/paraview_state_files/BST_2D.pvsm
:::
::::
:::::

:::::{tab-item} Full file
:sync: Full


```{literalinclude} ../../_static/gwb_input_files/BST_16_mass_conserving.wb
:language: json
:lineno-start: 1
:emphasize-lines: 37,39,40,42
```

::::{grid} 3
:::{grid-item-card} BST_16_mass_conserving.wb
:link: ../../_static/gwb_input_files/BST_16_mass_conserving.wb
:::
:::{grid-item-card} BST_16_mass_conserving.grid
:link: ../../_static/gwb_input_files/BST_16_mass_conserving.grid
:::
:::{grid-item-card} Paraview V3 state file 
:link: ../../_static/paraview_state_files/BST_v3.pvsm
:::
::::
:::::

::::::

## Changing to a mass conserving temperature model

One of the main ideas behind the mass conserving model is that the the slab not just warms up, but that the surrounding material cools down as well, thus conserving energy/mass. This means that the subducting plate feature needs to be able to change the temperature outside the slab. Given that a feature can only change the temperature inside it's defined area we need to expand that area. Hence, we need to increase the thickness which we originally set to 300 km. But this only extends the slab feature downwards, while we also want to cool down the material above the slab. For this we can use the `top truncation` parameter. The `top truncation` parameter is designed to make the slab thinner from above but when set to a negative value we can actually make it thicker! In this case we will set it to -100 km. 

Now that all the preparatory work is done, we can finally add the mass conserving model itself. The two most important new parameter are the `ridge coordinates` which should be the same as the subducting oceanic plate ridge coordinates in this case, and a parameter called `coupling depth`. The `coupling depth` defines the depth at which the slab surface first comes in contact with the hot mantle wedge.

::::::{tab-set}

:::::{tab-item} Important lines
:sync: Partial

```{literalinclude} ../../_static/gwb_input_files/BST_16_mass_conserving.wb
:language: json
:lineno-start: 44
:lines: 44-67
:emphasize-lines: 6,9,10,15,16,20,21,22,23
```
::::{grid} 3
:::{grid-item-card} BST_16_mass_conserving.wb
:link: ../../_static/gwb_input_files/BST_16_mass_conserving.wb
:::
:::{grid-item-card} BST_16_mass_conserving.grid
:link: ../../_static/gwb_input_files/BST_16_mass_conserving.grid
:::
:::{grid-item-card} Paraview 2D state file 
:link: ../../_static/paraview_state_files/BST_2D.pvsm
:::
::::
:::::

:::::{tab-item} Full file
:sync: Full


```{literalinclude} ../../_static/gwb_input_files/BST_16_mass_conserving.wb
:language: json
:lineno-start: 1
:emphasize-lines: 49,52,53,58,59,63,64,65,66
```

::::{grid} 3
:::{grid-item-card} BST_16_mass_conserving.wb
:link: ../../_static/gwb_input_files/BST_16_mass_conserving.wb
:::
:::{grid-item-card} BST_16_mass_conserving.grid
:link: ../../_static/gwb_input_files/BST_16_mass_conserving.grid
:::
:::{grid-item-card} Paraview V3 state file 
:link: ../../_static/paraview_state_files/BST_v3.pvsm
:::
::::
:::::

::::::



```{figure} ../../../../doc/sphinx/_static/images/user_manual/basic_starter_tutorial/BST_16.png
:name: BST_16
:alt: Basic Starter Tutorial section 18 highres result. 
:align: center

Basic Starter Tutorial section 16 high resolution result, where a mass conserving slab temperature is used. This has 8 times the resolution than the grid file above.
```

These are only the very basics of what the mass conserving temperature model can do. If you are interested in using this temperature model, please read the [Simple Subduction Model: 2D Cartesian](part:user_manual:chap:cookbooks:sec:simple_subduction_2d_cartesian)and [Simple Subduction Model: 2D Chunk](part:user_manual:chap:cookbooks:sec:simple_subduction_2d_chunk) cookbooks.
