
# The OPeNDAP extension to ASPECT

Provided here is the OPeNDAP client extension for ASPECT. 
OPeNDAP is both the name of a non-profit organization and the 
commonly-used name of a protocol which the OPeNDAP organization has 
developed. The OPeNDAP protocol provides a discipline-neutral means of 
requesting and providing data across the World Wide Web. OPeNDAP has
been working with the NSF/EarthCube funded BALTO (Brokered Alignment of 
Long-Tail Observations and Big Data) project to create a brokering 
cyberinfrastructure that facilitates access to unique data sources.
Using the BALTO/OPeNDAP extension of ASPECT will allow the user to query 
data from remote servers to be read in and used in local computations.

The OPeNDAP client extension to ASPECT has been tested on OSX (10.13,
10.14) and CentOS 7. It is currently available in source form only. In
this document you will find a description of the prerequisites for the
software, how to run the tests we have included, and how to make data
available from OPeNDAP servers so that ASPECT can use it.

While OPeNDAP makes its own open-source data server (called _Hyrax_),
that's not the only viable server that supports the _Data Access
Protocol_. Other servers include the THREDDS Data Server, ERDDAP, and
PyDAP. Any of these a can be used to provide data for use with ASPECT
(or any other software configured to work as a DAP client). For more
information, visit _www.opendap.org_, or the Unidata TDS, ERDDAP or
PyDAP sites (via Google).

## Building the software

To build the OPeNDAP extension to ASPECT, first get or build the
prerequisites and compile ASPECT with activated OPenDAP. The
prerequisites need to be installed before the ASPECT build is
configured so that the _cmake_ configuration tool will find them.

### Prerequisites

#### Install a binary copy of the libdap library.

On OSX, the easy way is to use homebrew:

> brew install libdap

The library is also available in Spack, see
https://spack.readthedocs.io/en/latest/package_list.html#libdap4

If you installed libdap using Spack or homebrew, skip down to
the next section.

#### Building libdap

If there are problems with any of the following steps, or you want to
install the software in a location other than the default, see the
detailed installation instructions. Skim over the instructions that
follow before you start to save some time.

The _libdap_ library has one prerequisite that is not commonly found
on OSX or Linux computers; it requires a recent version of the tool
_bison_. Run

> bison --version

and look at the version number. If it's less than 2.4 (or if there's
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

Get, build and install _libdap_ from source (but, of course, you don't
need to do this if you've installed a binary version of the library).

> wget https://www.opendap.org/pub/source/libdap-3.20.3.tar.gz
> tar -xzf libdap-3.20.3.tar.gz
> ./configure
> make
> sudo make install

### deal.II

To build ASPECT you will also need the deal.II library and it's
dependencies. as described in the ASPECT manual 
https://aspect-documentation.readthedocs.io/en/latest/index.html) and the wiki (https://github.com/geodynamics/aspect/wiki/).

### Compiling ASPECT

The only thing done differently to build the 'OPeNDAP enabled' version
of ASPECT is to provide information about _libdap_ to cmake when it
configures the ASPECT build using the ASPECT_WITH_LIBDAP option
(and optionally specify LIBDAP_DIR and LIBDAP_LIB) as follows:

> mkdir build
>
> cd build
>
> cmake -DLIBDAP_LIB=<libdap location>/build/lib/ -DLIBDAP_DIR=<libdap location>/build/ -DASPECT_WITH_LIBDAP=ON ..

Using -DASPECT_WITH_LIBDAP=OFF will make the client go back to running without
the libdap packaging (but not using libdap packages will make it impossible
for ASPECT to read data from a url). This option defaults to OFF meaning
libdap is optional.

## Running the tests

In the _build_ directory made above, run ASPECT using one of the
supplied _prm_ files. The file _opendap/prm\_files/aspect\_url\_test.prm_
will access data for lithospheric thickness for East Africa from a
data server running at OPeNDAP HQ; the _prm_ file named
_aspect\_test.prm_ (without _url_ in the name) will run the same
module using local data (those data can be found in
_opendap/input\_files/..._

> ./aspect -j ../opendap/prm\_files/aspect\_url\_test.prm

or

> mpirun -np 4 ./aspect ../opendap/prm\_files/aspect\_url\_test.prm

## Looking at the simulation

I've used _visit_ to render the simulation runs. It can be found here:
https://wci.llnl.gov/simulation/computer-codes/visit/.

