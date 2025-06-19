(sec:install:local-installation:fastscape)=
# Installing FastScape-fortran

To install FastScape-fortran, the first step is to add the following export to
the `.bashrc` file on linux or `.zprofile` on mac:

      export GFORTRAN_CONVERT_UNIT='big_endian'

This command tells FastScape how to write the VTK files. Without it the models will
run correctly, however FastScape visualizations will not show up in paraview. After this is
exported, next is to clone the code using the following command:

      git clone https://github.com/Djneu/fastscapelib-fortran

This branch includes modifications to the VTK writing for FastScape that are necessary for outputting within the ASPECT folders. You can download from the original repository as well, `https://github.com/fastscape-lem/fastscapelib-fortran`, however FastScape visualizations will be disabled in this case.

Next, run cmake with the option to build a shared library and that path to install fastscape to:

      cmake -DBUILD_FASTSCAPELIB_SHARED=ON \
             /path/to/fastscape_source

After this, compile the code using:

      make

Alternatively, you can pass a prefix using the command:

    cmake -DBUILD_FASTSCAPELIB_SHARED=ON \
                 -DCMAKE_INSTALL_PREFIX=/path/to/install/fastscape \
                 /path/to/fastscape_source

Then compile the code using:

      make install

This will create the cmake files in the current directory and install the required FastScape libraries inside `/path/to/install/fastscape/lib`.
