(part:user_manual:chap:cookbooks:sec:2d_cartesian_hydrated_slab)=
2D Hydrated Subducting Plate
======================

A simple setup of a subducting plate in a Cartesian box which contains hydrated lithologies. There are four lithologies, a 3 km
thick sediment layer, a 4 km thick mid-ocean ridge basalt (MORB) layer, a 4 thick gabbro layer, and a 9 km thick peridotite layer. The lithologies are specified to have an upper bound on the amount of bound water that they can hold of 3 wt%, 1 wt%, 0.5 wt%, and 2 wt% for sediment, MORB, gabbro, and peridotite, respectively. 

The relevant part of the World Builder file for prescribing the water content within the unsubducted oceanic plate looks like this:

:::{literalinclude} ../2d_cartesian_hydrated_slab.wb
:language: json
:lineno-start: 22
:lines: 22-36
:::

while the relevant part of the World Builder file for prescribing the water content within the subducting oceanic plate looks like this:

:::{literalinclude} ../2d_cartesian_hydrated_slab.wb
:language: json
:lineno-start: 38
:lines: 38-65
:::

And the bound water content within the subducting plate looks like this:

:::{figure} bound-water_isotherms.png
:name: slab_bound_water
:alt: 3D Cartesian Rift cookbook temperature contour. 
:align: center

The background is coloured by the bound wt% water, and isotherms are overlain on top at from 273 K to 1373 K at intervals of 100 K.
:::
