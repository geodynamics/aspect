
# Configuring ASPECT to Interface with External Software

For extended functionality, ASPECT can be configured with external
software. Below is a guide for setting the configuration to link
ASPECT with the external software that are currently supported.

# Geodynamic World Builder

The Geodynamic World Builder is a C++ tool for defining complex initial conditions for an ASPECT
simulation. By default, ASPECT will compile with the version of the Geodynamic World Builder
that comes bundled with the ASPECT source code, which is located in `$ASPECT_SOURCE_DIR/contrib/world_builder`.
However, to use a newer (or older) version of the Geodynamic World Builder, one can pass a flag to
`cmake` when configuring ASPECT. To do this, run the following command from the ASPECT build directory:

```
cmake -D WORLD_BUILDER_SOURCE_DIR=/path/to/world/builder/source ..
```

# FastScape

FastScape is a Fortran software for modeling surface processes. If you have downloaded and installed
FastScape-fortran (see {ref}`sec:install:local-installation:fastscape`), to use the coupling you must link ASPECT to the FastScape install.
This is done by running the following command from the ASPECT build directory:


    cmake -DASPECT_WITH_FASTSCAPE=ON \
                 -DFASTSCAPE_DIR=/path/to/fastscape/build ..

During cmake, a few outputs may be seen that notify if the link with FastScape is successful. First,
if the library is found the following will be seen:

     FastScape library found at /path/to/fastscape/build

If the FastScape version that allows visualizations with ASPECT is used, the following will be output:

     Found fastscape_named_vtk -- enabling its use

# Landlab

Landlab is a python library for modeling surface processes. This requires compiling ASPECT
with python. When compiling ASPECT with Landlab, we use the `uv` python package to ensure that a stable
python environment is always used, and that this environment always contains the necessary dependencies
for Landlab to run. For information on how to install `uv`, see the [uv documentation](https://docs.astral.sh/uv/).
Once uv is installed, to setup the correct python environment for landlab navigate to the $ASPECT_SOURCE_DIR
and run:

```
uv sync --project contrib/landlab/
```

Before compiling ASPECT with Landlab, activate this python environment by running:

```
source ./contrib/landlab/.venv/bin/activate
```

:::{warning}
Make sure that you deactivate any existing python/conda environments prior to compiling ASPECT with Landlab.
:::

With the environment configured, ASPECT can now be compiled against Landlab by running `cmake` with the following
flags from the ASPECT build directory:

```
cmake -D ASPECT_WITH_PYTHON=ON -D ASPECT_WITH_LANDLAB=ON -D Python3_EXECUTABLE=/path/to/ASPECT/SOURCE/DIR/contrib/landlab/.venv/bin/python3 ..
```
