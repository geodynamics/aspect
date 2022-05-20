
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

     cmake -DDEAL_II_DIR=/u/username/deal-installed/ ..

This will create an "out-of-source" build, where the build
directory is different from the source directory. While in-source builds
(where you run `cmake .` in your source directory) are supported, we strongly
recommend an out-of-source build as described above. Specifically, running the
whole test suite (see {ref}`sec:running_tests`) is only
supported this way.

:::{admonition} TODO
:class: error
Fix reference to section on running tests
:::
