PERPLEX INTEGRATION WITH ASPECT
===============================

About
-----

PerpleX is a code developed by Jamie Connolly which calculates
(amongst other things) equilibrium phase proportions and properties
given a bulk composition. It requires an endmember thermodynamic file,
a solution model file (several of which are bundled with PerpleX)
and a user-generated input file giving the names of these files,
the bulk composition of the rock, and the solutions and endmembers to
include/exclude during running of the various PerpleX programs.



Installation instructions
-------------------------

This directory should contain setup_perplex.bash and cmake_perplex.tar.
Run ./setup_perplex.sh, which should download the latest stable version
of the PerpleX source files, build them, and make a shared library for
use by ASPECT.



Preparing a PerpleX input file for ASPECT
-----------------------------------------

ASPECT uses a shared library which calls Gibbs minimization functions
inside PerpleX-meemum. In order to do this, a PerpleX input file is required
(an example named pyrolite.dat is given here).

If you know how to make a PerpleX input file already, great!
If not, PerpleX includes an executable called "build", which helps the user
build an input file step-by-step on the command-line.
You must know the relative path to the desired endmember
and solution model files.

ASPECT reassigns composition based on the values of compositional fields,
so it is largely unimportant which composition you choose when building
the input file. See the example parameter file named
perplex_lookup_composition.prm.



Troubleshooting / For more information
--------------------------------------

PerpleX has a website with detailed tutorials and information:
http://www.perplex.ethz.ch/, and an active discussion group:
http://www.perplex.ethz.ch/Perple_X_discussion_group.html.

If you get an error running one of the PerpleX executables, or you
don't know how to make a PerpleX input file, you should use the PerpleX
website/discussion group as a first port-of-call. If you get a PerpleX-related
error during your ASPECT run, please send an email to bob.myhill@bristol.ac.uk
or raise an issue on the ASPECT github page.



WARNING
-------

Please understand that the PerpleX Lookup material model is only a
proof-of-concept; the number of P-T-X evaluations is extremely large, which
means that ASPECT will be exceptionally slow for any real problems.
The development team is currently working on initialization routines which
will store PerpleX data before any calls to the material model evaluate
function. Stay tuned!
