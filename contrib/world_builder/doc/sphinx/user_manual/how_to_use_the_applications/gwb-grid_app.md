(part:user_manual:chap:how_to_use_the_apps:sec:gwb-grid_app)=
The gwb-grid app
================

While the `gwb-dat` app is great at providing individual points, providing gridded output can be very useful as well, especially when you want to to visualize the result. This is exactly what the `gwb-grid` app does. It creates a grid for you with the values you request and outputs it as a `vtu` file. The way of using `gwb-grid` is very similar to using `gwb-dat`: `gwb-grid  world_builder_file.wb grid_info.grid`. The only difference is that you now have to provide a `.grid` file instead of a `.dat` file.

A grid file consists of a number of parameters you can set. The start with a name (no spaces) followed by an equal sign, followed by the value. The available parameters are:

1. `grid_type`: Options are `cartesian`, `sphere`, `chunk` and `anullus`. 
2. `dim`: dimension, either 2 or 3.
3. `composition`: how many compositions
4. `vtu_output_format`: Either `ASCII` or readable output or `RawBinaryCompressed` for compressed output.
5. `x_min`: either the minimum x value or minimum longitude.
6. `x_max`: either the maximum x value or maximum longitude.
7. `y_min`: either the minimum y value or minimum longitude.
8. `y_min`: either the maximum y value or maximum longitude.
9. `z_min`: either the minimum z value or minimum radius.
10. `z_min`: either the maximum z value or maximum radius.
11. `n_cell_x`: either the cells in the x or longitude direction.
12. `n_cell_y`: either the cells in the y or latitude direction.
12. `n_cell_z`: either the cells in the z or radius direction.

An example of a grid file is the following: 

```{literalinclude} ../../../../tests/gwb-grid/spherical_subducting_plate_gridfile.grid 
:language: python
:lineno-start: 1
```

When you run it, it will produce a vtu file with the same name in the directory you run it from. 

More examples can be found in the `tests/gwb-grid/` directory.