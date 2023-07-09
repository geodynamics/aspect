
# Local installation

This is a brief explanation of how to compile and install the required
dependencies and ASPECT itself. This
installation procedure guarantees fastest runtimes, and largest flexibility,
but usually requires more work than the options mentioned in the previous
sections. While it is possible to install ASPECT's dependencies in
particular <span class="smallcaps">p4est</span>, <span
class="smallcaps">Trilinos</span>, and deal.II
manually, we recommend to use the `candi` software (see
<https://github.com/dealii/candi>). `candi` was written as an installation
program for deal.II, and includes a number of system specific instructions
that will be listed when starting the program. It can be flexibly configured
to allow for non-default compilers or libraries (e.g. Intel's MKL
instead of LAPACK) by changing entries in the configuration file `candi.cfg`,
or by providing platform specific installation files.

In case you encounter problems during the installation, please consult our
wiki (<https://github.com/geodynamics/aspect/wiki>) for frequently asked
questions and special instructions for MacOS users, before posting your
questions on the forum (<https://community.geodynamics.org/c/aspect>).


:::{toctree}
system-prereqs.md
using-candi.md
obtaining.md
compiling.md
documentation.md
:::
