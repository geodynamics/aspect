
# Obtaining ASPECT and initial configuration

The development version of ASPECT can be
downloaded by executing the command

     git clone https://github.com/geodynamics/aspect.git

If `$DEAL_II_DIR` points to your deal.II
installation, you can configure ASPECT by
running

     mkdir build; cd build; cmake ..

in the ASPECT directory created by the
`git clone` command above. If you did not set `$DEAL_II_DIR` you have to
supply cmake with the location:

     cmake -DDEAL_II_DIR=/home/username/deal-installed/ ..

This will create an "out-of-source" build, where the build
directory is different from the source directory, which is the
only supported option. (We do not support in-source builds where the source
and build directory are identical).

# Initial configuration with FastScape-fortran

If you have downloaded and installed FastScape-fortran, to use the coupling
you must link the ASPECT to the FastScape install. This is done by typing:


    cmake -DASPECT_WITH_FASTSCAPE=ON \
                 -DFASTSCAPE_DIR=/path/to/fastscape/build ..

During cmake, a few outputs may be seen that notify if the link with FastScape is successful. First,
if the library is found the following will be seen:

     FastScape library found at /path/to/fastscape/build

If the FastScape version that allows visualizations with ASPECT is used, the following will be output:

     Found fastscape_named_vtk -- enabling its use
