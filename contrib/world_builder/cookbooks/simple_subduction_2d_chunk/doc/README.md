(part:user_manual:chap:cookbooks:sec:simple_subduction_2d_chunk)=
Simple Subduction Model: 2D Chunk
=====================================
_Contributed by Magali Billen_

This tutorial will present how to set up a subduction model structure intended for use 
in a mantle convection model. The model will use a 2D slice in a spherical chunk geometry. 
The model will consist of 3 features:

1. An `oceanic plate` defining the overriding plate structure
2. An `oceanic plate` defining the sinking plate structure (along the model surface)
3. A `subducting plate` defining the subducted portion of the sinking plate (the slab)

For each feature, we need to choose the geometry and we can define the temperature model 
and the composition model to delineate layers in the plates and slab.

* All the files for this cookbook may be downloaded from the 
[github cookbook page](https://github.com/GeodynamicWorldBuilder/WorldBuilder/tree/main/cookbooks/simple_subduction_2d_chunk)
* It is maybe helpful to read the 
[Simple Subduction 2D Cartesian](part:user_manual:chap:cookbooks:sec:simple_subduction_2d_cartesian)
Cookbook before continuing with this cookbook.

For this example, the full model domain is a 90 degree profile in longitude and 1000 km 
deep.  The limits of the domain are chosen as 45 and 135 degrees, centered at 90 degrees. 
These limits are chosen in combination with the limits of the `cross section` parameter
in World Builder file (see below).  In the grid file, we choose a chunk geometry: in 3D, 
the chunk is a 3D section of a spherical shell. In 2D, the chunk is a great circle 
slice through a section of a spherical shell. We also indicate that there will be 2 compositions
in addition to the temperature.  

The figure below shows how the longitude angle is measured in the 3D World Builder space.  
the longitude parameter increases counter clockwise, from the y-axis. If viewed in Paraview,
this is the orientation in which the output will appear, . 

```{figure} chunk_2d_coordinate_sketch.png
:name: Chunk 2D geometry
:alt: Orientation of longitude coordinate in 2D chunk slice and location of the grid data points. 
:align: center
Orientation of longitude coordinate in 2D chunk slice and location of the grid data points. 
```

To visualize the output, using Paraview (or similar), we use two grid files
([gwb-grid](part:user_manual:chap:how_to_use_the_apps:sec:gwb-grid_app)):
a low resolution version (for quick checking) and a high resolution version (for careful 
checking of the details). The only parameters to change are the number of grid points
in each location:

:::::{tab-set}
::::{tab-item} Low resolution grid file 

:::{literalinclude} ../simple_chunk_2d_low_res.grid
:::
::::

::::{tab-item} High resolution grid file
:::{literalinclude} ../simple_chunk_2d_high_res.grid
:::

::::
:::::

These figures illustrate how the temperature structure may look _wrong_ (discontinuous) 
if the mesh resolution is too low, but increasing the mesh resolution, it is clear that 
the temperature is smooth across a fine mesh. 

::::{tab-set}

:::{tab-item} Slab with low resolution grid 
```{figure} chunk_2d_low_res_temperature.png
```
:::

:::{tab-item} Slab with high resolution grid 
```{figure} chunk_2d_high_res_temperature.png
```
:::

::::

 
## Root model information

In the root model section we need to indicate a spherical geometry (to agree with the
spherical coordinate system of the chunk). We also set the parameter 'depth method'
to "begin at end segment" 
(see [this page](part:user_manual:chap:concepts:sec:const_angle_spherical)) for more information
The `cross section` parameter defines the start and end locations for the cross section at 
the surface of the model (z = 0). For this example, we choose to put the cross section at 
the equator (y/latitude = 0), and give the start and end points of the cross section in 
longitude as x = 0 and x = 180 degrees. As noted above, the start value of the `cross section` 
can shift where features appear within the grid.

::::{tab-set}
:::{tab-item} Root model information
```{literalinclude} ../simple_subduction_2d_chunk.wb
:language: json
:lineno-start: 1
:lines: 1-11
```
:::
:::{tab-item} Full File
```{literalinclude} ../simple_subduction_2d_chunk.wb
:language: json
:lineno-start: 1
```
:::
::::

The three images below show how changing only the starting value for the cross
section shifts the position of the slab (no other values are changed; same grid file)
Therefore, unless you have a specific reason, the `cross section` parameter should 
start at x = 0 (longitude). For the models with a shifted start location, 
the temperature for the region beyond the sinking plate is filled in with the background 
adiabatic temperature gradient. 
 
:::::{tab-set}

::::{tab-item} Cross section: 0 to 180 
:::{figure} cross_section_0_180.png
:::
::::

::::{tab-item} Cross section: 10 to 180 
:::{figure} cross_section_10_180.png
:::
::::

::::{tab-item} Cross section: 30 to 180 
:::{figure} cross_section_30_180.png
:::
::::

:::::

The other physical parameters are defined as reference values, which will be used by the 
feature temperature models if the values are not entered separately.
 
### Feature Geometry, Temperature and Composition 

For each feature, we define the extent of that feature, a compositional layer and the
temperature structure.  For the compositional layer, we use a thickness of 100 km as
a marker of the outline of the plate, and not to indicate an actual compositional structure.
We also use two different compositional values: 0 for the overriding plate and 1 for the
sinking plate and subducting plate. 

#### Oceanic Plates
We will place the trench at the center of the profile at 90 degrees. Therefore,
* Overriding plate extends from 90 to 135 degrees in longitude.
* Sinking plate extends from 45 to 90 degrees in longitude. 

In the z direction, the maximum depth of the plate is chosen to be deep enough to account 
for any changes in the plate temperature in depth. If too small a value is chosen, a 
sharp jump in temperature will be seen in the output,
(see, for example, [Simple Subduction 2D Cartesian](part:user_manual:chap:cookbooks:sec:simple_subduction_2d_cartesian)).

For the temperature model we use the half space cooling model (appropriate for ages less
than 80 my). The temperature of the oceanic plates is calculated from the plate age, 
which is based on the location of the ridge for that plate and the spreading velocity. 
For this model we place the ridge for the sinking plate at the edge of the box (x = 135 degrees). 
However, for the overriding plate, we place the location of the ridge outside the model domain area 
(x = 0 degrees): this shows how World Builder considers distance (i.e., distance to 
the ridge) on a full 3D sphere, even though we are only interested in the effect within
the limits of the grid we define. 

To specify that the half space cooling model should meet the background adiabatic temperature, 
we set `bottom temperature` to -1. 
::::{tab-set}
:::{tab-item} Oceanic Plate information
```{literalinclude} ../simple_subduction_2d_chunk.wb
:language: json
:lineno-start: 12
:lines: 12-35
```
:::
:::{tab-item} Full File
```{literalinclude} ../simple_subduction_2d_chunk.wb
:language: json
:lineno-start: 1
```
:::
::::

#### Subducting Plate
Here we first add the compositional layer to continue the layer from the sinking plate. 
For the slab, we choose segments of different lengths, some with constant dip angle 
and others with variable dip. We also illustrate that the slab geometry can be overturned
by increasing the angle beyond 90 degree. We keep the `thickness` and `top truncation` 
values constant. These parameters are illustrated in 
[Simple Subduction 2D Cartesian](part:user_manual:chap:cookbooks:sec:simple_subduction_2d_cartesian)). 

The last step is to define the temperature structure of the slab. For this we choose the
temperature model `mass conserving` and the `half space model` option for `reference model name`.
The temperature structure depends on the age of the plate at the trench and subducting plate velocity. 
The age of the plate at the trench is determined from the `ridge coordinates`, and subducting 
`plate velocity`. The slab continues to get older with depth based on the length along the 
slab surface and subducting plate velocity. For the mass conserving temperature model, 
heating up of the slab and cooling down of the surrounding mantle depends on age, 
subducting plate velocity and length of the slab. 

In the temperature model `mass conserving`, the parameter `min distance slab top` 
has the same role as `top truncation` in defining how far above the slab layer 
information is calculated for the `subducting plate` feature. Similarly, in the temperature 
model `mass conserving`, the parameter `max distance slab top` has the same role 
as `thickness` in defining how far below the slab top information is calculated for 
the `subducting plate` feature. 

::::{tab-set}
:::{tab-item} Subducting Plate information
```{literalinclude} ../simple_subduction_2d_chunk.wb 
:language: json 
:lineno-start: 36
:lines: 36-62
```
:::
:::{tab-item} Full File
```{literalinclude} ../simple_subduction_2d_chunk.wb 
:language: json
:lineno-start: 1
```
:::
::::

Three other parameters control details of the temperature structure:
* the rate of heating is different where the slab is in contact with the overriding plate 
versus when it is in  direct contact with the mantle. This change is controlled by 
the `coupling depth` parameter.
* the temperature of the forearc region depends on the duration of subduction and
shielding of this region from the warmer mantle. The amount of additional cooling of the
forearc is controlled by the `forearc cooling factor`. This is a unconstrained tuning
parameter in this model. Typical values are in the range of 0-20.
* the down-dip end of the slab needs to smoothly transition into the surrounding mantle. 
To achieve this, the `taper distance`  indicates the distance at which to start linearly 
tapering the slab temperature into the background temperature.  

Final model figures:

:::::{tab-set} 

::::{tab-item} Full model domain 
:::{figure} chunk_2d_full_domain.png
Contours show outline of overriding plate composition (white) and sinking plate and slab
(yellow)
:::

::::

::::{tab-item} Center region of model domain
:::{figure} chunk_2d_high_res_temperature_compositions.png
Contours show outline of overriding plate composition (white) and sinking plate and slab
(yellow)
:::
::::

::::{tab-item} Corner of model 
:::{figure} chunk_2d_high_res_ridge_corner.png
Corner of model domain showing the ridge temperature structure for the sinking plate.
:::
::::

:::::

#### Experimenting with the Model
 
Recommendations for things to experiment with (only vary one parameter at a time):
* Try changing the shape of the slab, by adjusting the dips or lengths of the slab.
* Try change the location of one of the ridges.

 
