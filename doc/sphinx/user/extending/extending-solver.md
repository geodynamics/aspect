(sec:extending:solver)=
# Extending the basic solver

The core functionality of the code, i.e., that part of the code that
implements the time stepping, assembles matrices, solves linear and nonlinear
systems, etc., is in the `aspect::Simulator` class (see the [doxygen
documentation of this class](https://aspect.geodynamics.org/doc/doxygen/classaspect_1_1Simulator.html)).
Since the
implementation of this class has more than 3,000 lines of code, it is split
into several files that are all located in the `source/simulator` directory.
Specifically, functionality is split into the following files:

-   `source/simulator/core.cc`: This file contains the functions that drive
    the overall algorithm (in particular `Simulator::run`) through the main
    time stepping loop and the functions immediately called by
    `Simulator::run`.

-   `source/simulator/assembly.cc`: This is where all the functions are
    located that are related to assembling linear systems.

-   `source/simulator/solver.cc`: This file provides everything that has to do
    with solving and preconditioning the linear systems.

-   `source/simulator/initial_conditions.cc`: The functions in this file deal
    with setting initial conditions for all variables.

-   `source/simulator/checkpoint_restart.cc`: The location of functionality
    related to saving the current state of the program to a set of files and
    restoring it from these files again.

-   `source/simulator/helper_functions.cc`: This file contains a set of
    functions that do the odd thing in support of the rest of the simulator
    class.

-   `source/simulator/parameters.cc`: This is where we define and read
    run-time parameters that pertain to the top-level functionality of the
    program.

Obviously, if you want to extend this core functionality, it is useful to
first understand the numerical methods this class implements. To this end,
take a look at the paper that describes these methods, see
{cite:t}`kronbichler:etal:2012`. Further, there are two predecessor programs whose
extensive documentation is at a much higher level than the one typically found
inside ASPECT itself, since they are meant to
teach the basic components of convection simulators as part of the
deal.II tutorial:

-   The step-31 program at
    <https://www.dealii.org/developer/doxygen/deal.II/step_31.html>: This
    program is the first version of a convection solver. It does not run in
    parallel, but it introduces many of the concepts relating to the time
    discretization, the linear solvers, etc.

-   The step-32 program at
    <https://www.dealii.org/developer/doxygen/deal.II/step_32.html>: This is a
    parallel version of the step-31 program that already solves on a spherical
    shell geometry. The focus of the documentation in this program is on the
    techniques necessary to make the program run in parallel, as well as some
    of the consequences of making things run with realistic geometries,
    material models, etc.

Neither of these two programs is nearly as modular as
ASPECT, but that was also not the goal in creating
them. They will, however, serve as good introductions to the general approach
for solving thermal convection problems.

:::{note}
Neither this manual, nor the documentation in ASPECT makes much of an attempt at
teaching how to use the deal.II library upon which ASPECT is built. Nevertheless, you will
likely have to know at least the basics of deal.II to successfully work on the ASPECT code. We
refer to the resources listed at the beginning of this section as well as references {cite}`bangerth:etal:2007,bangerth:etal:2012`.
:::
