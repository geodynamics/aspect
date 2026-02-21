(part:GWB_parameter_listings:chap:visualization_grid_file:sec:index)=
Visualization Grid file
=======================

The visualization grid file is used to create a structured mesh and populate the node points with the generated world, i.e., initial conditions of temperature, composition, grains, etc. In essence, it is a text file that contains 3 major sections specifying the mesh points and values :

1. Mesh Output Type
* grid_type    = annulus|cartesian|chunk|spherical
The geometry of the output grid. The `gwb-grid` will not throw an error if the chosen grid type is inconsistent with the coordinate system in the input world_builder (`.wb`) file, but the generated output values would not make much sense. This means that for the spherical coordinate system, either of the `annulus`, `spherical`, or `chunk` grid types would give expected output for the computed properties.

* dim          = 2|3
The dimensions of the world. 

* compositions = 0
The number of compositions the output .vtu file contains. It can be any number from 0 to MAX_UNSIGNED_INT value. This means that even if we only have 1 composition in the input `.wb` file, it is possible to create visualization outputs for many compositional fields, which will be just equal to zero everywhere in the output mesh. If no value is specified, then it is assumed that there are no compositional fields, i.e., `compositions = 0`.

2. Mesh Extent
The following lines specify the domain of the grid. In case of 2D, you only need to specify the x and the z extents. Each of the following variables can be a double with a value from -MAX_DOUBLE to MAX_DOUBLE (inclusive). For the spherical coordinate system, the same variable names (as mentioned below) are used, where `x_{min|max}, y_{min|max}, z_{min|max}` represent the extent in the longitude, latitude, and radial directions, respectively.
* x_min
* x_max
* y_min
* y_max
* z_min 
* z_max

3. Mesh Discretization
Along with the domain extent, we also need to specify the number of points along each direction (for a 2D cross-section of world, only specify in the x and z directions). The following variables will do exactly that and can have any value from 1 to the MAX_INT (inclusive). For spherical coordinate system, the same variable names are used, where `n_cell_x, n_cell_y, n_cell_z` represent the discretization in the longitude, latitude, and radial directions, respectively.
* n_cell_x
* n_cell_y
* n_cell_z

An easier way to start writing your grid file is to copy any existing grid file available under the `tests` or `examples` directory and modify the file to include all the features from your `.wb` file.