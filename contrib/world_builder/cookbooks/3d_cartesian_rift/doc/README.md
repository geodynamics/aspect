(part:user_manual:chap:cookbooks:sec:3d_cartesian_rift)=
3D Cartesian Rift
======================

A simple setup of a rift in a Cartesian box, in which two plates are spreading at a constant spreading rate.
The files can be found [here](https://github.com/GeodynamicWorldBuilder/WorldBuilder/tree/main/cookbooks/3d_cartesian_rift).

The relevant part of the World Builder file looks like this:

:::{literalinclude} ../3d_cartesian_rift.wb
:language: json
:lineno-start: 1
:::

And the generated output model looks like this:

:::{figure} 1000K_temperature_contour.png
:name: 3D_cartesian_rift_temperature_contour
:alt: 3D Cartesian Rift cookbook temperature contour. 
:align: center

The 1000 K temperature isosurface of the 3D Cartesian rift cookbook. Colors represent the depth of the isosurface in meters and the
black box represents the boundaries of the model setup.
:::
