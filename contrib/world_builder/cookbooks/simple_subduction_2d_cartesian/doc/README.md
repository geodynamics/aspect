(part:user_manual:chap:cookbooks:sec:simple_subduction_2d_cartesian)=
Simple Subduction Model: 2D Cartesian
=====================================
_Contributed by Magali Billen_

This tutorial will present how to set up a subduction model structure intended for use 
in a mantle convection model. The model will use a 2D slice in a Cartesian geometry. 
The model will consist of 4 features:

1. A `mantle layer` with a uniform temperature
2. An `oceanic plate` defining the overriding plate structure
3. An `oceanic plate` defining the sinking plate structure (along the model surface)
4. A `subducting plate` defining the subducted portion of the sinking plate (the slab)

For each feature, we need to choose the geometry and we can define the temperature model 
and the composition model to delineate layers in the plates and slab.
 
For this tutorial, we will choose to use a simple temperature structure with a uniform 
background temperature (i.e., no adiabat) and we will use the `plate model` {cite}`Stein_Stein_1992` 
for temperature in the three plate related features. By using the same temperature model 
for all three features, it is possible to have continuous temperatures across the feature 
boundaries. 

When building a new World Builder model for input to some other software, we recommend 
creating a grid file so you can check the construction of your model as you add in 
additional features by visualizing the model output created by 
[gwb-grid](part:user_manual:chap:how_to_use_the_apps:sec:gwb-grid_app)
using Paraview (or similar). 
A grid file has a specific structure defining the geometry, extent of the model region 
to visualize and the grid spacing. Choose a grid spacing that is sufficient to check the 
model is defined properly. For a 2D model only the information for the 2D grid of points 
needs to be indicated in the grid file. For this example, the full model domain is 
8000 km long and 1600 km deep.  The grid file for this example is:

:::{literalinclude} ../simple_subduction_2d_cartesian.grid
:::

## Root model information

In the root model section we need to indicate the geometry (Cartesian). Even though we 
only want a 2D model, World Builder starts from a 3-D space. Therefore, we also need to 
define a 2-D cross section through the 3D World Builder space. The `cross section` 
parameter defines the start and end locations for the cross section at the surface of 
the model (z = 0). We choose to put the cross section at y = 0, and give the end 
points of the cross section at x = 0 and x = 8000e3 m. The other physical parameters 
are defined as reference values, which will be used by the feature temperature models 
if the values are not entered separately.

::::{tab-set}
:::{tab-item} Root model information
```{literalinclude} ../simple_subduction_2d_cartesian.wb
:language: json
:lineno-start: 1
:lines: 1-6
```
:::
:::{tab-item} Full File
```{literalinclude} ../simple_subduction_2d_cartesian.wb
:language: json
:lineno-start: 1
```
:::
::::
 
### Feature Geometry, temperature and composition 

Because the World Builder defines features in 3D space, even though we want to define 
a 2D model, we still need to indicate some lateral extent for features in the y-direction. 
Any value can be chosen for this width, so we will center the cross section at y = 0, and 
let the y-direction extend from -100e3 m to 100e3 m.

#### Mantle Layer
We start by defining the mantle layer to encompass all of the model domain. 
We don't need to define a thickness for the mantle, but we need to define its lateral 
extent. This is done by assigning the x, y coordinates of the vertices of a polygon 
given in clockwise order. In this example, the polygon is a rectangle, so there are 4 
pairs of coordinates provide.  Because this is the first feature in the parameter file, 
features that follow and overlap with this layer will replace (by default) the existing value. 

The reason we define a mantle layer is so we can set the background temperature to a u
uniform value of 1573 K. The same value will be used as the bottom temperature in the 
plate models. A constant temperature mantle is a simplification, but is sometimes useful 
for testing ideas. A large model domain size is chosen to reduce the effects of the side 
boundaries, However, an even wider model is likely needed with a depth of only 1600 km. 
Alternatively, the box depth can be increased to 3000 km.  

::::{tab-set}
:::{tab-item} Root model information
```{literalinclude} ../simple_subduction_2d_cartesian.wb
:language: json
:lineno-start: 7
:lines: 7-11
```
:::
:::{tab-item} Full File
```{literalinclude} ../simple_subduction_2d_cartesian.wb
:language: json
:lineno-start: 1
```
:::
::::

#### Oceanic Plates
The geometry of the oceanic plates is defined by their lateral extent in the x direction
and their depth in the z direction, again setting the y direction width to an arbitrary value: 
* Overriding plate extends from 0 to 3500e3 m in the x-direction.
* Sinking plate extends from 3500e3 m to 8000e3 m in the x-direction. 

In the z direction, the maximum depth of the plate is chosen to be deep enough to account 
for any changes in the plate temperature in depth. If too small a value is chosen, a 
sharp jump in temperature will be seen in the output. The down side to using a larger value
is that it can take longer to needlessly assign values to additional grid points in the 
larger depth range. After defining the temperature model, by trial and error, we choose the
maximum depth to be 200e3 m for the overriding plate and 300e3 m for the sinking plate. 

::::{tab-set}
:::{tab-item} max depth 200 km
```{figure} oceanic_plates_max_depth_too_small.png
:name: max depth too small
:alt: Parameter max depth = 200 km is too small for subducting plate 
:align: center
Parameter max depth = 200 km is too small for subducting plate 
```
:::
:::{tab-item} max depth 300 km
```{figure} oceanic_plates_max_depth.png
:name: max depth 300 km
:alt: Parameter max depth = 300 km is large enough for age of subducting plate
:align: center
Parameter max depth = 300 km is large enough for age of subducting plate
```
:::
::::

For the temperature model we use the plate model: for plate ages less than about 80 my the
half space {cite}`turcotte_schubert_2014` and plate models  {cite}`Stein_Stein_1992` 
are similar, but for older ages, the half space model over predicts the thickness of the 
lithosphere, based on bathymetric data.  The temperature of the oceanic plates is calculated 
from the plate age, which is based on the location of the ridge for that plate and 
the spreading velocity. For this model we place the ridge for each plate at the box edges 
at x = 0 and x = 8000e3 m. 

Finally, a composition field has also been included for the sinking plate. A thickness
of 50e3 m is chosen so this layer can easily be seen in the example. This layer does 
not correspond to a geological feature, but can be used to refer to or show the initial 
location or thickness of the slab.  The same approach can be used to define a crustal 
layer with a max depth/thickness of 8e3 m. 

::::{tab-set}
:::{tab-item} Root model information
```{literalinclude} ../simple_subduction_2d_cartesian.wb
:language: json
:lineno-start: 12
:lines: 12-31
```
:::
:::{tab-item} Full File
```{literalinclude} ../simple_subduction_2d_cartesian.wb
:language: json
:lineno-start: 1
```
:::
::::

#### Subducting Plate
Here we illustrate how to define a generic slab profile by defining segments of the slab
extending sequentially deeper from the surface trench location, given by the `coordinates`.
The shape of the slab is controlled by the length of each segment and the dip angle of the 
segment, which can vary linearly along the length of each segment. Each segment also has
a `thickness`, which has the same purpose as the max depth parameter for the oceanic plate 
feature, and a `top truncation` parameter, which sets the distance above the slab surface
over which the composition or temperature will also be calculated. Distances above the slab
are negative. The code block below illustrates how the parameters defining the segments are
input in the World Builder file. 

::::{tab-set}
:::{tab-item} Root model information
```{literalinclude} ../simple_subduction_2d_cartesian.wb
:language: json
:lineno-start: 32
:lines: 32-41
```
:::
:::{tab-item} Full File
```{literalinclude} ../simple_subduction_2d_cartesian.wb
:language: json
:lineno-start: 1
```
:::
::::

Here we first add the compositional layer for visualization (above) of the location and 
thickness of the slab. We chose segments of different lengths, some with constant dip angle 
and others with variable dip. We keep the `thickness` and `top truncation` values constant. 
The figure below also visually illustrates the parameters defining the segments. 

```{figure} segment_geometry_composition.png
:name: slab segment geometry
:alt: Composition layer for the slab illustrating the slab segment geometry information.
:align: center
Composition layer for the slab illustrating the slab segment geometry information.
```

The last step is to define the temperature structure of the slab. For this we choose the
temperature model `mass conserving` and the `plate model` option for `reference model name`.
Similar to the plate model temperature structure for the oceanic plate, the temperature
structure depends on the age of the plate at the trench and subducting plate velocity. 
The age of the plate at the trench is determined from the `ridge coordinates`, and subducting 
`plate velocity`. The slab continues to get older with depth based on the length along the 
slab surface and subducting plate velocity. For the mass conserving temperature model, 
heating up of the slab and cooling down of the surrounding mantle depends on age, 
subducting plate velocity and length of the slab. The parameter `min distance slab top` 
has the same role as `top truncation` and the parameter `max distance slab top` has the 
same role as `thickness`.

::::{tab-set}
:::{tab-item} Root model information

```{literalinclude} ../simple_subduction_2d_cartesian.wb
:language: json
:lineno-start: 42
:lines: 42-51
```
:::
:::{tab-item} Full File
```{literalinclude} ../simple_subduction_2d_cartesian.wb
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

::::{tab-set}
:::{tab-item} Final model full domain
```{figure} simple_subduction_2d_cartesian_full.png
:name: Full model
:alt: Final temperature for 2d cartesian simple subduction model, full model domain
:align: center
Final temperature for 2d cartesian simple subduction model, full model domain
```
:::
:::{tab-item} Center of model domain
```{figure} simple_subduction_2d_cartesian_center.png
:name: Center of model
:alt: Final temperature for 2d cartesian simple subduction model, center of model
:align: center
Final temperature for 2d cartesian simple subduction model, center of model
```
:::
::::

#### Experimenting with the Model
 
Recommendations for things to experiment with (only vary one parameter at a time) - change the:
* `max depth` for the oceanic plate: make smaller until its too small 
* `top truncation`: make smaller until its too small
* `coupling depth`: make deeper or shallower - how does temperature change 
* `forearc cooling factor`: make larger or smaller - how does temperature change
* `taper distance`: make larger or smaller - how does slab tip change
* `spreading velocity` for the overriding plate: make smaller - how does temperature change
* `spreading velocity` for the sinking plate: make smaller - how does temperature change 
* `plate velocity` for the slab to match that of the sinking plate
 
