
# System prerequisites

`candi` will show system specific instructions on startup, but its
prerequisites are relatively widely used and packaged for most operating
systems. You will need compilers for C, C++ and Fortran, the GNU make system,
the CMake build system, and the libraries and header files of BLAS, LAPACK and
zlib, which is used for compressing the output data. To use more than one
process for your computations, you will need to install a MPI library, its
headers and the necessary executables to run MPI programs. There are some
optional packages for additional features, like the HDF5 libraries for
additional output formats and Numdiff for checking
ASPECT's test results with reasonable accuracy,
but these are not strictly required, and in some operating systems, they are
not available as packages but need to be compiled from scratch. Finally, for
obtaining a recent development version of
ASPECT, you will need the git version control system.

An exemplary command to obtain all required packages on Ubuntu 20.04 would be:

    sudo apt-get install build-essential \
                         cmake \
                         gcc \
                         g++ \
                         gfortran \
                         git \
                         libblas-dev \
                         liblapack-dev \
                         libopenmpi-dev \
                         numdiff \
                         openmpi-bin \
                         zlib1g-dev
