# Using Aspect with Fastscapelib (C++)

## Install Fastscapelib using conda

You can use conda to install fastscapelib and all of its dependencies (we assume
that you are familiar with managing conda environments):

```
conda install fastscapelib -c conda-forge
```

## Install Fastscapelib from source

You can follow the steps below to install Fastscapelib and its dependencies from
source.

### Set the installation prefix path for Fastscapelib and its dependencies

Choose a common location where you want to download the source files of
Fastscapelib and its dependencies:

```
$ export FS_SOURCE_PATH=/path/to/fastscapelib/source
```

Choose a common location where you want to install Fastscapelib and its
dependencies:

```
export CMAKE_INSTALL_PREFIX=/path/to/fastscapelib/install/dir
```

### Install xtl and xtensor

Download the `xtl` git repository and checkout the last stable version:

```
$ cd $FS_SOURCE_PATH
$ git clone https://github.com/xtensor-stack/xtl
$ cd xtl
$ git checkout 0.7.7 
```

Build and install `xtl`:

```
$ cmake -S. -Bbuild
$ cmake --install build
```

Download the `xtensor` git repository and checkout the last stable version:

```
$ cd $FS_SOURCE_PATH
$ git clone https://github.com/xtensor-stack/xtensor
$ cd xtensor
$ git checkout 0.25.0 
```

Build and install `xtensor`:

```
$ cmake -S. -Bbuild
$ cmake --install build
```

### Install Fastscapelib

Download the `fastscapelib` git repository and checkout the last stable version:

```
$ cd $FS_SOURCE_PATH
$ git clone https://github.com/fastscape-lem/fastscapelib
$ cd fastscapelib
$ git checkout v0.2.2
```

IMPORTANT: edit the `cmake/fastscapelibConfig.cmake` file with the editor of
your choice and add the following line just below `@PACKAGE_INIT@` (this will be
required every time you checkout another version or branch and re-install
Fastscapelib, but this won't be required for Fastscapelib version >0.2.2):

```
include(CMakeFindDependencyMacro)
```

Build and install `fastscapelib`:

```
$ cmake -S. -Bbuild
$ cmake --install build
```

## Compile Aspect with Fastscapelib enabled

If you have installed Fastscapelib using conda, don't forget to activate the
environment where you have installed Fastscapelib. Once activated, configuring
`aspect` with CMake should then find Fastscapelib automatically.

If you have installed Fastscapelib from source like described above, you can
help CMake find the Fastscapelib installation path using:

```
export CMAKE_PREFIX_PATH=/path/to/fastscapelib/install/dir
```

where the given path corresponds to the path that you set during installation
via the `CMAKE_INSTALL_PREFIX` environment variable (see instructions above).

Configuring `aspect` should then find Fastscapelib. Check the output of the
following command run from aspect's source root directory (note: you might need
to set additional arguments or run extra commands in order to find the other
dependencies of `aspect`, we skipped it here):

```
cmake -S. -Bbuild/temp -DASPECT_WITH_FASTSCAPELIB=ON
```
