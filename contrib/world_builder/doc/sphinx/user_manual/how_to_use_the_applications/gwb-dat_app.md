(part:user_manual:chap:how_to_use_the_apps:sec:gwb-dat_app)=
The gwb-dat app
=====================

The `gwb-dat` app allows you to query the GWB world at many points at a time. As arguments, it takes a World Builder file, which specifies what the world looks like, and a `.dat` file, which specifies what points you want to query and what information you want. An example of how it is run is the following: `gwb-dat world_builder_file.wb point_info.dat`. There is only one optional extra argument which is an option to limit some consistency checks in debug mode: `--limit-debug-consistency-checks` (or `-ldcc` as a shortcut). If you want to use it, it must be added after the `.dat` file.

The main part of the `.dat` file consists of 3 (in 2D) or 4 (in 3D) columns of numbers. The first 2 or 3 columns are the coordinates of the 2D or 3D point respectively. The last column is the depth. In debug mode, spherical models are checked by default for  consistency. You can turn this off using the above mentioned `--limit-debug-consistency-checks` option. 

Lines with a `#` are either ignored or can have a special meaning. For example, a line with `# random text here` is ignored, but a line with `# dim = 2` sets the dimension to 2. The available options are:

1. `dim = `: 2 or 3
2. `compositions = `: number of compositions
3. `grain compositions = `: number of grain compositions
4. `number of grains = `: number of grains per particle.
5. `convert spherical = `: true or false. This option is only allowed to be set true in 3D. It allows you to input radius, longitude and latitude instead of x, y and z.

An example of a `.dat` file is the following:

```{literalinclude} ../../../../tests/gwb-dat/app_oceanic_plate_cartesian.dat
:language: python
:lineno-start: 1
```

This file can be found in the tests directory: `tests/gwb-dat/app_oceanic_plate_cartesian.dat`.

Running `gwb-dat tests/gwb-dat/app_oceanic_plate_cartesian.wb tests/gwb-dat/app_oceanic_plate_cartesian.dat` will give the following output: 

```{literalinclude} ../../../../tests/gwb-dat/app_oceanic_plate_cartesian/screen-output.log
:language: python
:lineno-start: 1
```

You can see it appends the temperature and 9 compositions. It also provides a header line with a symbol indicating what each field means.

More examples of `.dat` files can be found in the `tests/gwb-dat/` directory.
