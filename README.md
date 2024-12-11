ASPECT - the Advanced Solver for Planetary Evolution, Convection, and Tectonics
====================================================================================

[![License GPL2+:](https://img.shields.io/badge/License-GPL%202%2B-red)](https://github.com/geodynamics/aspect/blob/main/LICENSE)
[![Online Documentation](https://readthedocs.org/projects/aspect-documentation/badge/?version=latest)](https://aspect-documentation.readthedocs.io/en/latest/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14371679.svg)](https://doi.org/10.5281/zenodo.14371679)
[![pdf manual](https://img.shields.io/badge/get-PDF-green.svg)](https://doi.org/10.6084/m9.figshare.4865333)
[![Docker image:](https://github.com/geodynamics/aspect/actions/workflows/docker.yml/badge.svg)](https://hub.docker.com/r/geodynamics/aspect)
[![Tester image:](https://github.com/geodynamics/aspect/actions/workflows/build_dealii_master.yml/badge.svg)](https://github.com/geodynamics/aspect/actions/workflows/build_dealii_master.yml)

About
-----

ASPECT is a code to simulate convection and tectonic processes in
Earth and other planetary bodies. It has grown from a pure
mantle-convection code into a tool for many geodynamic applications
including applications for inner core convection, lithospheric scale
deformation, two-phase flow, and numerical methods development. The
project is supported by CIG (https://geodynamics.org).


Installation instructions
-------------------------

The steps to install the necessary dependencies and ASPECT itself are described
in the Installation instructions section of the ASPECT
[manual](https://aspect-documentation.readthedocs.io/en/latest/user/install/index.html). If you encounter
problems during the installation, please consult our
[wiki](https://github.com/geodynamics/aspect/wiki) for typical installation
problems or specific instructions for MacOS users, before asking your question
on the forum.

In short, ASPECT is configured using CMake and has the following requirements:
- CMake 3.13.4
- GCC, Clang, or Intel compiler with C++14 support
- [deal.II](https://github.com/dealii/dealii) 9.5.0 or newer configured with:
  - MPI, Trilinos, p4est, SUNDIALS (required)
  - BLAS/LAPACK, zlib (strongly recommended)
  - HDF5 (optional)
- optional: PerpleX, LIBDAP, NETCDF

Running and extending ASPECT
----------------------------

Instructions on how to run and extend, as well as on how to interpret the
output of ASPECT can also be found in the ASPECT
[manual](https://aspect-documentation.readthedocs.io/en/latest/index.html). This manual also
discusses the structure of the source code.

For getting started, you can also watch recordings from our
[tutorials](https://github.com/geodynamics/aspect/wiki/Tutorial-Materials).


Contributing to ASPECT
----------------------

ASPECT is a community project that lives by the participation of its
members — i.e., including you! It is our goal to build an inclusive
and participatory community so we are happy that you are interested in
participating! We have collected a set of guidelines and advice on how
to get involved in the community and keep them in the
[CONTRIBUTING.md](CONTRIBUTING.md)
file in ASPECT's repository.



More information
----------------

For more information see:
 - The official website at https://aspect.geodynamics.org
 - The current [Online Documentation](https://aspect-documentation.readthedocs.io/en/latest/)
 - [How to cite ASPECT](https://aspect.geodynamics.org/cite.html)
 - For questions on the source code of ASPECT, portability, installation, new or existing features, etc., use the [ASPECT forum](https://community.geodynamics.org/c/aspect). This forum is where the ASPECT users and developers all hang out.
 - See the [Future Plan Board](https://github.com/geodynamics/aspect/projects/2) for a high level overview of items in development.
 - ASPECT is primarily based on the deal.II library. If you have particular questions about deal.II, contact the [deal.II discussion groups](https://www.dealii.org/mail.html).
 - ASPECT is being developed by a large, collaborative, and inclusive community. It is currently maintained by the following people:
     - Wolfgang Bangerth: bangerth@math.colostate.edu
     - Juliane Dannberg: judannberg@gmail.com
     - Menno Fraters: mrfraters@ucdavis.edu
     - Rene Gassmoeller: rene.gassmoeller@mailbox.org
     - Anne Glerum: acglerum@gfz-potsdam.de
     - Timo Heister: heister@clemson.edu
     - Bob Myhill: bob.myhill@bristol.ac.uk
     - John Naliboff: john.naliboff@nmt.edu

 - The following people have significantly contributed and furthered ASPECT's goals and are therefore Principal Developers:

   - Jacky Austermann
   - Wolfgang Bangerth
   - Juliane Dannberg
   - Menno Fraters
   - Rene Gassmoeller
   - Anne Glerum
   - Timo Heister
   - Bob Myhill
   - John Naliboff
   - Cedric Thieulot

 - A complete and growing list of the many authors that have contributed over the years can be found at [GitHub contributors](https://github.com/geodynamics/aspect/graphs/contributors).
 - If you have specific questions about ASPECT that are not suitable for public and archived mailing lists, feel free to contact the maintainers or principal developers.



License
-------

ASPECT is published under [GPL v2 or newer](LICENSE).
