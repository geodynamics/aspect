
# The OPeNDAP extension to ASPECT

The OPeNDAP client extension to ASPECT has been tested on OSX
(10.13, 10.14) and is currently available in source form only. In this
document you will find a description of the prerequisites for the
software, how to run the tests we have included and how to make data
available from OPeNDAP servers so that ASPECT can use it.

While OPeNDAP makes its own open-source data server (called _Hyrax_),
that's certainly not the only viable server that supports the _Data
Access Protocol_. Other servers include the THREDDS Data Server,
ERDDAP, and PyDAP. Any of these a can be used to provide data for use
with ASPECT (or any other software configured to work as a DAP
client). For more information, visit _www.opendap.org_, or the Unidata
TDS, ERDDAP or PyDAP sites (via Google).

## Building the software

To build the OPeNDAP extension to ASPECT, first build the prerequisites
and then then the modified version of ASPECT. The prerequisites need to
be built before the ASPECT build is configured the  _cmake_
configuration tool will find them.

### Prerequisites - Building libdap

If there are problems with any of the following steps, or you want to
install the software in a location other than the default, see the
detailed installation instructions. SKim over the instructions that
follow before you start to save some time.

The _libdap_ library has one prerequisite that is not commonly found
on OSX or Linux computers; it requires a recent version of the tool
_bison_. Run

> bison --version

and look at the version number. If it's less than 3.1 (or if there's
no such program installed on your computer), download, unpack, build
and install a recent version of _bison_ from
https://ftp.gnu.org/gnu/bison/. We have used version 3.2 and
version 3.3 These instructions will accomplish this:

> wget https://ftp.gnu.org/gnu/bison/bison-3.3.tar.xz
> tar -xzf bison-3.3.tar.xz
> cd bison-3.3
> ./configure
> make
> sudo make install

Get and install _libdap_

> wget https://www.opendap.org/pub/source/libdap-3.20.3.tar.gz
> tar -xzf libdap-3.20.3.tar.gz
> ./configure
> make
> sudo make install

### DEAL.II

To build ASPECT you will also need the DEAL.II library and it's
dependencies.

#### OSX
For OSX, go to https://www.dealii.org/ and follow the Download
button/link to get a pre-built dmg for OSX. To install the binary,
open the dmg and drag the 'deal.II' icon to the Applications
directory.

#### CentOS7
For CentOS 7, Build using candi as per the ASPECT manual. For me that meant 
getting candi from GitHub (git clone https://github.com/aspect/candi)
and running

> ./candi.sh --platform=deal.II-toolchain/platforms/supported/centos7.platform

and then following its directions before proceeding with the actual
build. That meant installing a list of packages:

> sudo yum install patch svn git wget @development-tools gcc-c++ cmake \
> openmpi openmpi-devel patch  libtool libtool-ltdl libtool-ltdl-devel \
> lua lua-devel blas blas-devel lapack lapack-devel  doxygen graphviz \
> graphviz-devel qt-devel

Then get a new shell and run:

> module load mpi/openmpi-x86_64

(apparently the 'module' tcl program is part of @development-tools. You
need to get a new shell for the command to work.)

Then you need to add the following env vars to the shell:

> export CC=mpicc; export CXX=mpicxx; export FC=mpif90; export FF=mpif77

Now actually run candi (use "./candi --platform=deal.II-toolchain/platforms
/supported/centos7.platform" again) this time hitting return when prompted
to continue (assuming you have done the various things). You will need to
run candi once, make a list, do the stuff and then run it again.
  
### Compiling ASPECT

The only thing done differently to build the 'OPeNDAP enabled' version
of ASPECT is to provide information about _libdap_ to cmake when it
configures the ASPECT build using the LIBDAP\_DIR and
LIBDAP\_ON options combined with the DEAL\_II\_DIR option as follows:

> mkdir build
> cd build
> cmake .. -DLIBDAP\_DIR=/usr/local/ -DLIBDAP=ON
>               -DDEAL\_II\_DIR=/Applications/deal.II-9.0.0.app/Contents/Resources/

## Running the tests

In the _build_ directory made above, run ASPECT using one of the
supplied _prm_ files. The file _opendap/prm\_files/aspect\_url\_test.prm_
will access data for lithospheric thickness for East Africa from a
data server running at OPeNDAP HQ; the _prm_ file named
_aspect\_test.prm_ (without _url_ in the name) will run the same
modle using local data (those data can be found in
_opendap/input\_files/..._

> ./aspect -j ../opendap/prm\_files/aspect\_url\_test.prm

or

> mpirun -np 4 ./aspect ../opendap/prm\_files/aspect\_url\_test.prm

## Looking at the simulation

I've used _visit_ to render the simulation runs. It can be found here:
https://wci.llnl.gov/simulation/computer-codes/visit/.

## Making new data accessible to ASPECT using an OPeNDAP server

## About the modifications

## Detailed information about building
