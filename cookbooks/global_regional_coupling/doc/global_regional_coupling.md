(sec:cookbooks:global_regional_coupling)=
# Prescribing Velocities from Global Convection Models as Boundary Conditions in Regional Models

*This section was contributed by Daniel Douglas.*

In this cookbook we showcase how to use the python script `extract_local_velocity.py` (located in the [/aspect/contrib/python/] directory) to prescribe velocity boundary conditions on a regional spherical chunk using the velocity solution from a global convection model. For the `extract_local_velocity.py` script to work, pvpython must be installed. The idea is that when simulating a 3D spherical chunk the choice of boundary conditions has a non-negligible effect on the regional solution. One can negate the boundary effects by increasing the domain of the regional model, but the far-field effects from global tectonics are still generally ignored by the regional model. By prescribing the velocity solution from a global convection model, the boundary conditions import some of the context of far-field global tectonics that would otherwise be missing from the regional model. We will compare the difference in the resulting solution from an instantaneous Stokes solution in a regional chunk when applying velocity boundary conditions from a global convection model, or applying free slip boundary conditions. In both cases, we will use the S20RTS model as the initial temperature distribution, and the simple material model. We choose the spherical chunk domain spanning from -55$^\circ$ < latitude < -20$\circ$, 152$\circ$ < longitude < 210$\circ$, and 4770 km < radius < 6370 km.

To start, output from a global convection model must be generated. We will run the `initial-condition-S20RTS` cookbook for this demonstration, since this is a relatively low resolution and doesn't require much computational effort. Assuming ASPECT has been compiled, we can run this cookbook using the following command:

```
$ASPECT_SOURCE_DIR/./aspect-release $ASPECT_SOURCE_DIR/cookbooks/initial-conditions-S20RTS/S20RTS.prm
```

This will generate a `solution.pvd` file, so we can now run the `extract_local_velocity.py` script by first navigating to the correct directory:

```
cd $ASPECT_SOURCE_DIR/contrib/python
pvpython extract_local_velocity.py
```

The script will write four ASCII files, one for each of the lateral boundaries (North, East, South, and West) into a directory called `regional_velocity_files`, which are configured so that they can be directly used in an ASPECT model. The `extract_local_velocity.py` script can be modified to output along different slices by changing the `latitude_bounds`, `longitude_bounds`, and `radius_bounds` parameters. For this case, plotting the ASCII files quickly in python shows that the global velocity along each boundary looks like this:

```{figure-md} fig:pvpython-output
<img src="pvpython_output.png" style="width:90.0%" />

The velocities extracted from the output of the `initial-conditions-S20RTS` global convection model.
```

 Now that the velocity boundary conditions are ready to be used in ASPECT, we can compare how the exact same regional model is impacted by the inclusion of far-field effects outside of the boundaries of the regional model. The input file `free_slip_boundaries.prm` assumes all six boundaries are closed, which is achieved by defining them as free slip boundaries:

```{literalinclude} free_slip_boundaries.part.prm
```

whereas the input file `global_regional_coupling.prm` applies the ASCII files written by `extract_local_velcoity.py` on the North, East, South and West boundary, while applying an open bottom boundary and a free slip top boundary.

```{literalinclude} global_regional_coupling.part.prm
```

Specifying an open bottom boundary is important in these models because it helps the model maintain mass conservation. While mass is certainly conserved within a global convection model, any given subset of that global model is not required to have a flow field that conserves mass. Allowing one model boundary to be open gives the regional model room to adjust it's flow pattern to ensure that mass is conserved within the subset of the global model. We compare the velocities along each of the 6 boundaries, as well as the internal flow field, in the series of figures below.

```{figure-md} fig:top-boundary
<img src="top_boundary.png" style="width:90.0%" />

Comparison of the velocity on the top boundary of the regional model with free slip boundary conditions (left) and with global velocity boundary conditions (right).
```

```{figure-md} fig:bottom-boundary
<img src="bottom_boundary.png" style="width:90.0%" />

Comparison of the velocity on the bottom boundary of the regional model with free slip boundary conditions (left) and with global velocity boundary conditions (right).
```

```{figure-md} fig:east-boundary
<img src="east_boundary.png" style="width:90.0%" />

Comparison of the velocity on the East boundary of the regional model with free slip boundary conditions (left) and with global velocity boundary conditions (right).
```

```{figure-md} fig:west-boundary
<img src="west_boundary.png" style="width:90.0%" />

Comparison of the velocity on the West boundary of the regional model with free slip boundary conditions (left) and with global velocity boundary conditions (right).
```

```{figure-md} fig:south-boundary
<img src="south_boundary.png" style="width:90.0%" />

Comparison of the velocity on the South boundary of the regional model with free slip boundary conditions (left) and with global velocity boundary conditions (right).
```

```{figure-md} fig:north-boundary
<img src="north_boundary.png" style="width:90.0%" />

Comparison of the velocity on the North boundary of the regional model with free slip boundary conditions (left) and with global velocity boundary conditions (right).
```

Visually, it is clear that the velocity magnitudes are greatly impacted by the application of the global velocity boundary conditions. Also, the lateral boundaries for the model which applies the global velocity boundary condition recovers the expected velocity which we plotted in Figure {numref}`fig:pvpython-output`. It is also visually clear that the internal flow patterns themselves are dramatically different, showcasing how much of an impact the inclusion of far field effects can be on the direction and magnitude of flow within regional models.
