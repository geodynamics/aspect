(sec:install:local-installation:compiling)=
# Compiling ASPECT

After downloading ASPECT and having built the
libraries it builds on, you can compile it by typing

      make

on the command line (or `make -jN` if you have multiple processors in your
machine, where `N` is the number of processors). This builds the
ASPECT executable which will reside in the `build`
directory and will be named `aspect` or `aspect-release`.
To run ASPECT from the main source directory, you would need
to reference it as `./build/aspect`.
