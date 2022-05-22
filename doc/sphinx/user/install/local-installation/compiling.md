(sec:install:local-installation:compiling)=
# Compiling ASPECT and generating documentation

After downloading ASPECT and having built the
libraries it builds on, you can compile it by typing

      make

on the command line (or `make -jN` if you have multiple processors in your
machine, where `N` is the number of processors). This builds the
ASPECT executable which will reside in the `build`
directory and will be named `aspect`. To run
ASPECT from the main source directory, you would need
to reference it as `./build/aspect`. If you intend to modify
ASPECT for your own experiments, you may want to also
generate documentation about the source code. This can be done using the
command

      cd doc; make

which assumes that you have the `doxygen` documentation generation tool
installed. Most Linux distributions have packages for `doxygen`. The result
will be the file `doc/doxygen/index.html` that is the starting point for
exploring the documentation.
