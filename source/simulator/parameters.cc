/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/simulator.h>
#include <aspect/global.h>
#include <aspect/utilities.h>
#include <aspect/melt.h>
#include <aspect/volume_of_fluid/handler.h>
#include <aspect/newton.h>
#include <aspect/mesh_deformation/free_surface.h>

#include <deal.II/base/parameter_handler.h>

#include <sys/stat.h>
#include <stdlib.h>
#include <boost/lexical_cast.hpp>

namespace aspect
{
  template <int dim>
  Parameters<dim>::Parameters (ParameterHandler &prm,
                               MPI_Comm mpi_communicator)
  {
    parse_parameters (prm, mpi_communicator);
  }


  template <int dim>
  void
  Parameters<dim>::
  declare_parameters (ParameterHandler &prm)
  {
    prm.declare_entry ("Dimension", "2",
                       Patterns::Integer (2,3),
                       "The number of space dimensions you want to run this program in. "
                       "ASPECT can run in 2 and 3 space dimensions.");
    prm.declare_entry ("Additional shared libraries", "",
                       Patterns::List (Patterns::FileName()),
                       "A list of names of additional shared libraries that should be loaded "
                       "upon starting up the program. The names of these files can contain absolute "
                       "or relative paths (relative to the directory in which you call ASPECT). "
                       "In fact, file names that do not contain any directory "
                       "information (i.e., only the name of a file such as <myplugin.so> "
                       "will not be found if they are not located in one of the directories "
                       "listed in the \\texttt{LD_LIBRARY_PATH} environment variable. In order "
                       "to load a library in the current directory, use <./myplugin.so> "
                       "instead."
                       "\n\n"
                       "The typical use of this parameter is so that you can implement "
                       "additional plugins in your own directories, rather than in the ASPECT "
                       "source directories. You can then simply compile these plugins into a "
                       "shared library without having to re-compile all of ASPECT. See the "
                       "section of the manual discussing writing extensions for more "
                       "information on how to compile additional files into a shared "
                       "library.");

    prm.declare_entry ("Resume computation", "false",
                       Patterns::Selection ("true|false|auto"),
                       "A flag indicating whether the computation should be resumed from "
                       "a previously saved state (if true) or start from scratch (if false). "
                       "If auto is selected, models will be resumed if there is an existing "
                       "checkpoint file, otherwise started from scratch.");

    prm.declare_entry ("Max nonlinear iterations", "10",
                       Patterns::Integer (1),
                       "The maximal number of nonlinear iterations to be performed.");

    prm.declare_entry ("Max nonlinear iterations in pre-refinement", boost::lexical_cast<std::string>(std::numeric_limits<int>::max()),
                       Patterns::Integer (0),
                       "The maximal number of nonlinear iterations to be performed in the pre-refinement "
                       "steps. This does not include the last refinement step before moving to timestep 1. "
                       "When this parameter has a larger value than max nonlinear iterations, the latter is used.");

    prm.declare_entry ("Start time", "0.",
                       Patterns::Double (),
                       "The start time of the simulation. Units: Years if the "
                       "'Use years in output instead of seconds' parameter is set; "
                       "seconds otherwise.");

    prm.declare_entry ("Timing output frequency", "100",
                       Patterns::Integer(0),
                       "How frequently in timesteps to output timing information. This is "
                       "generally adjusted only for debugging and timing purposes. If the "
                       "value is set to zero it will also output timing information at the "
                       "initiation timesteps.");

    prm.declare_entry ("Use years in output instead of seconds", "true",
                       Patterns::Bool (),
                       "When computing results for mantle convection simulations, "
                       "it is often difficult to judge the order of magnitude of results "
                       "when they are stated in MKS units involving seconds. Rather, "
                       "some kinds of results such as velocities are often stated in "
                       "terms of meters per year (or, sometimes, centimeters per year). "
                       "On the other hand, for non-dimensional computations, one wants "
                       "results in their natural unit system as used inside the code. "
                       "If this flag is set to `true' conversion to years happens; if "
                       "it is `false', no such conversion happens. Note that when `true', "
                       "some input such as prescribed velocities should also use years "
                       "instead of seconds.");

    prm.declare_entry ("CFL number", "1.0",
                       Patterns::Double (0.),
                       "In computations, the time step $k$ is chosen according to "
                       "$k = c \\min_K \\frac {h_K} {\\|u\\|_{\\infty,K} p_T}$ where $h_K$ is the "
                       "diameter of cell $K$, and the denominator is the maximal magnitude "
                       "of the velocity on cell $K$ times the polynomial degree $p_T$ of the "
                       "temperature discretization. The dimensionless constant $c$ is called the "
                       "CFL number in this program. For time discretizations that have explicit "
                       "components, $c$ must be less than a constant that depends on the "
                       "details of the time discretization and that is no larger than one. "
                       "On the other hand, for implicit discretizations such as the one chosen "
                       "here, one can choose the time step as large as one wants (in particular, "
                       "one can choose $c>1$) though a CFL number significantly larger than "
                       "one will yield rather diffusive solutions. Units: None.");

    prm.declare_entry ("Maximum time step",
                       /* boost::lexical_cast<std::string>(std::numeric_limits<double>::max() /
                                                           year_in_seconds) = */ "5.69e+300",
                       Patterns::Double (0.),
                       "Set a maximum time step size for the solver to use. Generally the time step "
                       "based on the CFL number should be sufficient, but for complicated models "
                       "or benchmarking it may be useful to limit the time step to some value. "
                       "The default value is a value so that when converted from years into seconds "
                       "it equals the largest number representable by a floating "
                       "point number, implying an unlimited time step."
                       "Units: Years or seconds, depending on the ``Use years "
                       "in output instead of seconds'' parameter.");

    prm.declare_entry ("Maximum first time step",
                       /* boost::lexical_cast<std::string>(std::numeric_limits<double>::max() /
                                                           year_in_seconds) = */ "5.69e+300",
                       Patterns::Double (0.),
                       "Set a maximum time step size for only the first timestep. Generally the time step "
                       "based on the CFL number should be sufficient, but for complicated models "
                       "or benchmarking it may be useful to limit the first time step to some value, "
                       "especially when using the free surface, which needs to settle to prevent "
                       "instabilities. This should in that case be combined with a value set for "
                       "``Maximum relative increase in time step''. "
                       "The default value is a value so that when converted from years into seconds "
                       "it equals the largest number representable by a floating "
                       "point number, implying an unlimited time step. "
                       "Units: Years or seconds, depending on the ``Use years "
                       "in output instead of seconds'' parameter.");

    prm.declare_entry ("Maximum relative increase in time step", boost::lexical_cast<std::string>(std::numeric_limits<int>::max()),
                       Patterns::Double (0.),
                       "Set a percentage with which the the time step is limited to increase. Generally the "
                       "time step based on the CFL number should be sufficient, but for complicated models "
                       "which may suddenly drastically change behavior, it may be useful to limit the increase "
                       "in the time step, without limiting the time step size of the whole simulation to a "
                       "particular number. For example, if this parameter is set to $50$, then that means that "
                       "the time step can at most increase by 50\\% from one time step to the next, or by a "
                       "factor of 1.5. "
                       "Units: \\%.");

    prm.declare_entry ("Use conduction timestep", "false",
                       Patterns::Bool (),
                       "Mantle convection simulations are often focused on convection "
                       "dominated systems. However, these codes can also be used to "
                       "investigate systems where heat conduction plays a dominant role. "
                       "This parameter indicates whether the simulator should also use "
                       "heat conduction in determining the length of each time step.");

    const std::string allowed_solver_schemes = "single Advection, single Stokes|iterated Advection and Stokes|"
                                               "single Advection, iterated Stokes|no Advection, iterated Stokes|"
                                               "no Advection, single Stokes|"
                                               "no Advection, iterated defect correction Stokes|"
                                               "single Advection, iterated defect correction Stokes|"
                                               "iterated Advection and defect correction Stokes|"
                                               "iterated Advection and Newton Stokes|single Advection, iterated Newton Stokes|"
                                               "single Advection, no Stokes|IMPES|iterated IMPES|"
                                               "iterated Stokes|Newton Stokes|Stokes only|Advection only|"
                                               "first timestep only, single Stokes|no Advection, no Stokes";

    prm.declare_entry ("Nonlinear solver scheme", "single Advection, single Stokes",
                       Patterns::Selection (allowed_solver_schemes),
                       "The kind of scheme used to resolve the nonlinearity in the system. "
                       "`single Advection, single Stokes' means that no nonlinear iterations are done, "
                       "and the temperature, compositional fields and Stokes equations are solved exactly "
                       "once per time step, one after the other. "
                       "The `iterated Advection and Stokes' scheme iterates this decoupled approach "
                       "by alternating the solution of the temperature, composition and Stokes systems. "
                       "The `single Advection, iterated Stokes' scheme solves the temperature and composition "
                       "equation once at the beginning of each time step and then iterates out the solution of "
                       "the Stokes equation. "
                       "The `no Advection, iterated Stokes' scheme only solves the Stokes system, iterating "
                       "out the solution, and ignores compositions and the temperature equation (careful, "
                       "the material model must not depend on the temperature or composition; this is mostly "
                       "useful for Stokes benchmarks). "
                       " The `no Advection, single Stokes' scheme only solves the Stokes system once per "
                       "timestep. This is also mostly useful for Stokes benchmarks. "
                       "The `single Advection, no Stokes' scheme only solves the temperature and other advection "
                       "systems once, and instead of solving for the Stokes system, a prescribed velocity "
                       "and pressure is used. "
                       "The `iterated Advection and Newton Stokes' scheme iterates by alternating the solution "
                       "of the temperature, composition and Stokes equations, using Picard iterations for the "
                       "temperature and composition, and Newton iterations for the Stokes system. "
                       "The `single Advection, iterated Newton Stokes' scheme solves "
                       "the temperature and composition equations once at the beginning of each time step and "
                       "then iterates out the solution of the Stokes equation, using Newton iterations for the "
                       "Stokes system. "
                       "The `first timestep only, single Stokes' scheme solves the Stokes equations exactly "
                       "once, at the first time step. No nonlinear iterations are done, and the temperature and "
                       "composition systems are not solved. "
                       "\n\n"
                       "The `IMPES' scheme is deprecated and only allowed for reasons of backwards "
                       "compatibility. It is the same as `single Advection, single Stokes' ."
                       "The `iterated IMPES' scheme is deprecated and only allowed for reasons of "
                       "backwards compatibility. It is the same as `iterated Advection and Stokes'. "
                       "The `iterated Stokes' scheme is deprecated and only allowed for reasons of "
                       "backwards compatibility. It is the same as `single Advection, iterated Stokes'. "
                       "The `Stokes only' scheme is deprecated and only allowed for reasons of "
                       "backwards compatibility. It is the same as `no Advection, iterated Stokes'. "
                       "The `Advection only' scheme is deprecated and only allowed for reasons of "
                       "backwards compatibility. It is the same as `single Advection, no Stokes'. "
                       "The `Newton Stokes' scheme is deprecated and only allowed for reasons of "
                       "backwards compatibility. It is the same as `iterated Advection and Newton Stokes'.");

    prm.declare_entry ("Nonlinear solver tolerance", "1e-5",
                       Patterns::Double(0., 1.),
                       "A relative tolerance up to which the nonlinear solver will iterate. "
                       "This parameter is only relevant if the `Nonlinear solver scheme' does nonlinear "
                       "iterations, in other words, if it is set to something other than "
                       "`single Advection, single Stokes' or `single Advection, no Stokes'.");

    prm.declare_entry ("Pressure normalization", "surface",
                       Patterns::Selection ("surface|volume|no"),
                       "If and how to normalize the pressure after the solution step. "
                       "This is necessary because depending on boundary conditions, "
                       "in many cases the pressure is only determined by the model "
                       "up to a constant. On the other hand, we often would like to "
                       "have a well-determined pressure, for example for "
                       "table lookups of material properties in models "
                       "or for comparing solutions. If the given value is `surface', then "
                       "normalization at the end of each time steps adds a constant value "
                       "to the pressure in such a way that the average pressure at the surface "
                       "of the domain is what is set in the `Surface pressure' parameter; "
                       "the surface of the domain is determined by asking "
                       "the geometry model whether a particular face of the geometry has a zero "
                       "or small `depth'. If the value of this parameter is `volume' then the "
                       "pressure is normalized so that the domain average is zero. If `no' is "
                       "given, the no pressure normalization is performed.");

    prm.declare_entry ("Surface pressure", "0.",
                       Patterns::Double(),
                       "The value the pressure is normalized to in each time step when "
                       "`Pressure normalization' is set to `surface' with default value 0. "
                       "This setting is ignored in all other cases."
                       "\n\n"
                       "The mathematical equations that describe thermal convection "
                       "only determine the pressure up to an arbitrary constant. On "
                       "the other hand, for comparison and for looking up material "
                       "parameters it is important that the pressure be normalized "
                       "somehow. We do this by enforcing a particular average pressure "
                       "value at the surface of the domain, where the geometry model "
                       "determines where the surface is. This parameter describes what "
                       "this average surface pressure value is supposed to be. By "
                       "default, it is set to zero, but one may want to choose a "
                       "different value for example for simulating only the volume "
                       "of the mantle below the lithosphere, in which case the surface "
                       "pressure should be the lithostatic pressure at the bottom "
                       "of the lithosphere."
                       "\n\n"
                       "For more information, see the section in the manual that discusses "
                       "the general mathematical model.");

    prm.declare_entry ("Adiabatic surface temperature", "0.",
                       Patterns::Double(),
                       "In order to make the problem in the first time step easier to "
                       "solve, we need a reasonable guess for the temperature and pressure. "
                       "To obtain it, we use an adiabatic pressure and temperature field. "
                       "This parameter describes what the `adiabatic' temperature would "
                       "be at the surface of the domain (i.e. at depth zero). Note "
                       "that this value need not coincide with the boundary condition "
                       "posed at this point. Rather, the boundary condition may differ "
                       "significantly from the adiabatic value, and then typically "
                       "induce a thermal boundary layer."
                       "\n\n"
                       "For more information, see the section in the manual that discusses "
                       "the general mathematical model.");

    prm.declare_entry ("Output directory", "output",
                       Patterns::DirectoryName(),
                       "The name of the directory into which all output files should be "
                       "placed. This may be an absolute or a relative path.");

    prm.declare_entry ("Use operator splitting", "false",
                       Patterns::Bool(),
                       "If set to true, the advection and reactions of compositional fields and "
                       "temperature are solved separately, and can use different time steps. Note that "
                       "this will only work if the material/heating model fills the reaction\\_rates/"
                       "heating\\_reaction\\_rates structures. Operator splitting can be used with any "
                       "existing solver schemes that solve the temperature/composition equations.");

    prm.declare_entry ("World builder file", "",
                       Patterns::FileName(),
                       "Name of the world builder file. If empty, the world builder is not initialized.");

    prm.enter_subsection ("Solver parameters");
    {
      prm.declare_entry ("Temperature solver tolerance", "1e-12",
                         Patterns::Double(0., 1.),
                         "The relative tolerance up to which the linear system for "
                         "the temperature system gets solved. See `Stokes solver "
                         "parameters/Linear solver tolerance' for more details.");

      prm.declare_entry ("Composition solver tolerance", "1e-12",
                         Patterns::Double(0., 1.),
                         "The relative tolerance up to which the linear system for "
                         "the composition system gets solved. See `Stokes solver "
                         "parameters/Linear solver tolerance' for more details.");

      prm.enter_subsection ("Advection solver parameters");
      {
        prm.declare_entry ("GMRES solver restart length", "50",
                           Patterns::Integer(1),
                           "This is the number of iterations that define the "
                           "GMRES solver restart length. Increasing this "
                           "parameter makes the solver more robust and decreases "
                           "the number of iterations. Be aware that "
                           "increasing this number increases the memory usage "
                           "of the advection solver, and makes individual "
                           "iterations more expensive.");
      }
      prm.leave_subsection();


      prm.enter_subsection ("Stokes solver parameters");
      {
        prm.declare_entry ("Stokes solver type", "block AMG",
                           Patterns::Selection(StokesSolverType::pattern()),
                           "This is the type of solver used on the Stokes system. The block geometric "
                           "multigrid solver currently has a limited implementation and therefore "
                           "may trigger Asserts in the code when used. If this is the case, "
                           "please switch to 'block AMG'. Additionally, the block GMG solver requires "
                           "using material model averaging.");

        prm.declare_entry ("Use direct solver for Stokes system", "false",
                           Patterns::Bool(),
                           "If set to true the linear system for the Stokes equation will "
                           "be solved using Trilinos klu, otherwise an iterative Schur "
                           "complement solver is used. The direct solver is only efficient "
                           "for small problems.");

        prm.declare_entry ("Krylov method for cheap solver steps", "GMRES",
                           Patterns::Selection(StokesKrylovType::pattern()),
                           "This is the Krylov method used to solve the Stokes system. Both options, GMRES "
                           "and IDR(s), solve non-symmetric, indefinite systems. GMRES "
                           "guarantees the residual will be reduced in each iteration while IDR(s) has "
                           "no such property. On the other hand, the vector storage requirement for GMRES is "
                           "dependent on the restart length and can be quite restrictive (since, for "
                           "the matrix-free GMG solver, memory is dominated by these vectors) whereas "
                           "IDR(s) has a short term recurrence. "
                           "Note that the IDR(s) Krylov method is not available for the AMG solver since "
                           "it is not a flexible method, i.e., it cannot handle a preconditioner which "
                           "may change in each iteration (the AMG-based preconditioner contains a CG solve "
                           "in the pressure space which may have different number of iterations each step).");

        prm.declare_entry ("IDR(s) parameter", "2",
                           Patterns::Integer(1),
                           "This is the sole parameter for the IDR(s) Krylov solver and will dictate the "
                           "number of matrix-vector products and preconditioner applications per iteration (s+1) "
                           "and the total number of temporary vectors required (5+3*s). For s=1, this method is "
                           "analogous to BiCGStab. As s is increased this method is expected to converge to "
                           "GMRES in terms of matrix-vector/preconditioner applications to solution.");

        prm.declare_entry ("Linear solver tolerance", "1e-7",
                           Patterns::Double(0., 1.),
                           "A relative tolerance up to which the linear Stokes systems in each "
                           "time or nonlinear step should be solved. The absolute tolerance will "
                           "then be $\\| M x_0 - F \\| \\cdot \\text{tol}$, where $x_0 = (0,p_0)$ "
                           "is the initial guess of the pressure, $M$ is the system matrix, "
                           "$F$ is the right-hand side, and tol is the parameter specified here. "
                           "We include the initial guess of the pressure "
                           "to remove the dependency of the tolerance on the static pressure. "
                           "A given tolerance value of 1 would "
                           "mean that a zero solution vector is an acceptable solution "
                           "since in that case the norm of the residual of the linear "
                           "system equals the norm of the right hand side. A given "
                           "tolerance of 0 would mean that the linear system has to be "
                           "solved exactly, since this is the only way to obtain "
                           "a zero residual."
                           "\n\n"
                           "In practice, you should choose the value of this parameter "
                           "to be so that if you make it smaller the results of your "
                           "simulation do not change any more (qualitatively) whereas "
                           "if you make it larger, they do. For most cases, the default "
                           "value should be sufficient. In fact, a tolerance of 1e-4 "
                           "might be accurate enough.");

        prm.declare_entry ("Number of cheap Stokes solver steps", "200",
                           Patterns::Integer(0),
                           "As explained in the paper that describes ASPECT (Kronbichler, Heister, and Bangerth, "
                           "2012, see \\cite{KHB12}) we first try to solve the Stokes system in every "
                           "time step using a GMRES iteration with a poor but cheap "
                           "preconditioner. By default, we try whether we can converge the GMRES "
                           "solver in 200 such iterations before deciding that we need a better "
                           "preconditioner. This is sufficient for simple problems with variable "
                           "viscosity and we never need the second phase with the more expensive "
                           "preconditioner. On the other hand, for more complex problems, and in "
                           "particular for problems with strongly nonlinear viscosity, the 200 "
                           "cheap iterations don't actually do very much good and one might skip "
                           "this part right away. In that case, this parameter can be set to "
                           "zero, i.e., we immediately start with the better but more expensive "
                           "preconditioner.");

        prm.declare_entry ("Maximum number of expensive Stokes solver steps", "1000",
                           Patterns::Integer(0),
                           "This sets the maximum number of iterations used in the expensive Stokes solver. "
                           "If this value is set too low for the size of the problem, the Stokes solver will "
                           "not converge and return an error message pointing out that the user didn't allow "
                           "a sufficiently large number of iterations for the iterative solver to converge.");

        prm.declare_entry ("GMRES solver restart length", "50",
                           Patterns::Integer(1),
                           "This is the number of iterations that define the GMRES solver restart length. "
                           "Increasing this parameter helps with convergence issues arising from high localized "
                           "viscosity jumps in the domain. Be aware that increasing this number increases the "
                           "memory usage of the Stokes solver, and makes individual Stokes iterations more "
                           "expensive.");

        prm.declare_entry ("Linear solver A block tolerance", "1e-2",
                           Patterns::Double(0., 1.),
                           "A relative tolerance up to which the approximate inverse of the $A$ block "
                           "of the Stokes system is computed. This approximate $A$ is used in the "
                           "preconditioning used in the GMRES solver. The exact definition of this "
                           "block preconditioner for the Stokes equation can be found in "
                           "\\cite{KHB12}.");

        prm.declare_entry ("Use full A block as preconditioner", "false",
                           Patterns::Bool(),
                           "This parameter determines whether we use an simplified approximation of the "
                           "$A$ block as preconditioner for the Stokes solver, or the full $A$ block. The "
                           "simplified approximation only contains the terms that describe the coupling of "
                           "identical components (plus boundary conditions) as described in "
                           "\\cite{KHB12}. The full block is closer to the description in "
                           "\\cite{rudi2017weighted}."
                           "\n\n"
                           "There is no clear way to determine which preconditioner "
                           "performs better. The default value (simplified approximation) requires "
                           "more outer GMRES iterations, but is faster to apply in each iteration. The full block "
                           "needs less assembly time (because the block is "
                           "available anyway), converges in less GMRES iterations, but requires more time per "
                           "iteration. There are also differences in the amount of memory consumption between "
                           "the two approaches."
                           "\n\n"
                           "The default value should be good for relatively simple models, but in "
                           "particular for very strong viscosity contrasts the full $A$ block can be "
                           "advantageous.");

        prm.declare_entry ("Linear solver S block tolerance", "1e-6",
                           Patterns::Double(0., 1.),
                           "A relative tolerance up to which the approximate inverse of the $S$ block "
                           "(i.e., the Schur complement matrix $S = BA^{-1}B^{T}$) of the Stokes "
                           "system is computed. This approximate inverse of the $S$ block is used "
                           "in the preconditioning used in the GMRES solver. The exact definition of "
                           "this block preconditioner for the Stokes equation can be found in "
                           "\\cite{KHB12}.");
      }
      prm.leave_subsection ();

      prm.enter_subsection ("AMG parameters");
      {
        prm.declare_entry ("AMG smoother type", "Chebyshev",
                           Patterns::Selection ("Chebyshev|symmetric Gauss-Seidel"),
                           "This parameter sets the type of smoother for the AMG. "
                           "The default is strongly recommended for any normal runs "
                           "with ASPECT. There are some indications that the symmetric "
                           "Gauss-Seidel might be better and more stable for the Newton "
                           "solver. For extensive benchmarking of various settings of the "
                           "AMG parameters in this section for the Stokes problem and others, "
                           "see https://github.com/geodynamics/aspect/pull/234.");

        prm.declare_entry ("AMG smoother sweeps", "2",
                           Patterns::Integer(0),
                           "Determines how many sweeps of the smoother should be performed. When the flag elliptic "
                           "is set to true, (which is true for ASPECT), the polynomial degree of "
                           "the Chebyshev smoother is set to this value. The term sweeps refers to the number of "
                           "matrix-vector products performed in the Chebyshev case. In the non-elliptic case, "
                           "this parameter sets the number of SSOR relaxation sweeps for post-smoothing to be performed. "
                           "The default is strongly recommended. There are indications that for the Newton solver a different "
                           "value might be better. For extensive benchmarking of various settings of the "
                           "AMG parameters in this section for the Stokes problem and others, "
                           "see https://github.com/geodynamics/aspect/pull/234.");

        prm.declare_entry ("AMG aggregation threshold", "0.001",
                           Patterns::Double(0., 1.),
                           "This threshold tells the AMG setup how the coarsening should be performed. "
                           "In the AMG used by ML, all points that strongly couple with the tentative coarse-level "
                           "point form one aggregate. The term strong coupling is controlled by the variable "
                           "aggregation\\_threshold, meaning that all elements that are not smaller than "
                           "aggregation\\_threshold times the diagonal element do couple strongly. "
                           "The default is strongly recommended. There are indications that for the Newton solver a different "
                           "value might be better. For extensive benchmarking of various settings of the "
                           "AMG parameters in this section for the Stokes problem and others, "
                           "see https://github.com/geodynamics/aspect/pull/234.");

        prm.declare_entry ("AMG output details", "false",
                           Patterns::Bool(),
                           "Turns on extra information on the AMG solver. Note that this will generate much more output.");
      }
      prm.leave_subsection ();
      prm.enter_subsection ("Operator splitting parameters");
      {
        prm.declare_entry ("Reaction time step", "1000.0",
                           Patterns::Double (0.),
                           "Set a time step size for computing reactions of compositional fields and the "
                           "temperature field in case operator splitting is used. This is only used "
                           "when the nonlinear solver scheme ``operator splitting'' is selected. "
                           "The reaction time step must be greater than 0. "
                           "If you want to prescribe the reaction time step only as a relative value "
                           "compared to the advection time step as opposed to as an absolute value, you "
                           "should use the parameter ``Reaction time steps per advection step'' and set "
                           "this parameter to the same (or larger) value as the ``Maximum time step'' "
                           "(which is 5.69e+300 by default). "
                           "Units: Years or seconds, depending on the ``Use years "
                           "in output instead of seconds'' parameter.");

        prm.declare_entry ("Reaction time steps per advection step", "0",
                           Patterns::Integer (0),
                           "The number of reaction time steps done within one advection time step "
                           "in case operator splitting is used. This is only used if the nonlinear "
                           "solver scheme ``operator splitting'' is selected. If set to zero, this "
                           "parameter is ignored. Otherwise, the reaction time step size is chosen according to "
                           "this criterion and the ``Reaction time step'', whichever yields the "
                           "smaller time step. "
                           "Units: none.");
      }
      prm.leave_subsection ();
      prm.enter_subsection ("Diffusion solver parameters");
      {
        prm.declare_entry ("Diffusion length scale", "1.e4",
                           Patterns::Double (0.),
                           "Set a length scale for the diffusion of compositional fields if the "
                           "``prescribed field with diffusion'' method is selected for a field. "
                           "More precisely, this length scale represents the square root of the "
                           "product of diffusivity and time in the diffusion equation, and controls "
                           "the distance over which features are diffused. "
                           "Units: \\si{\\meter}.");
      }
      prm.leave_subsection ();
    }
    prm.leave_subsection ();

    prm.enter_subsection("Formulation");
    {
      prm.declare_entry ("Formulation", "custom",
                         Patterns::Selection ("isentropic compression|custom|anelastic liquid approximation|Boussinesq approximation"),
                         "Select a formulation for the basic equations. Different "
                         "published formulations are available in ASPECT (see the list of "
                         "possible values for this parameter in the manual for available options). "
                         "Two ASPECT specific options are\n"
                         "\\begin{enumerate}\n"
                         "  \\item `isentropic compression': ASPECT's original "
                         "formulation, using the explicit compressible mass equation, "
                         "and the full density for the temperature equation.\n"
                         "  \\item `custom': A custom selection of `Mass conservation' and "
                         "`Temperature equation'.\n"
                         "\\end{enumerate}\n\n"
                         "\\note{Warning: The `custom' option is "
                         "implemented for advanced users that want full control over the "
                         "equations solved. It is possible to choose inconsistent formulations "
                         "and no error checking is performed on the consistency of the resulting "
                         "equations.}\n\n"
                         "\\note{The `anelastic liquid approximation' option here can also be "
                         "used to set up the `truncated anelastic liquid approximation' as long as "
                         "this option is chosen together with a material model that defines a "
                         "density that depends on temperature and depth and not on the pressure.}");

      prm.declare_entry ("Mass conservation", "ask material model",
                         Patterns::Selection ("incompressible|isentropic compression|hydrostatic compression|"
                                              "reference density profile|implicit reference density profile|"
                                              "projected density field|"
                                              "ask material model"),
                         "Possible approximations for the density derivatives in the mass "
                         "conservation equation. Note that this parameter is only evaluated "
                         "if `Formulation' is set to `custom'. Other formulations ignore "
                         "the value of this parameter.");
      prm.declare_entry ("Temperature equation", "real density",
                         Patterns::Selection ("real density|reference density profile"),
                         "Possible approximations for the density in the temperature equation. "
                         "Possible approximations are `real density' and `reference density profile'. "
                         "Note that this parameter is only evaluated "
                         "if `Formulation' is set to `custom'. Other formulations ignore "
                         "the value of this parameter.");
      prm.declare_entry ("Enable additional Stokes RHS", "false",
                         Patterns::Bool (),
                         "Whether to ask the material model for additional terms for the right-hand side "
                         "of the Stokes equation. This feature is likely only used when implementing force "
                         "vectors for manufactured solution problems and requires filling additional outputs "
                         "of type AdditionalMaterialOutputsStokesRHS.");
      prm.declare_entry ("Enable elasticity", "false",
                         Patterns::Bool (),
                         "Whether to include the additional elastic terms on the right-hand side of "
                         "the Stokes equation.");
      prm.declare_entry ("Enable prescribed dilation", "false",
                         Patterns::Bool (),
                         "Whether to include additional terms on the right-hand side of "
                         "the Stokes equation to set a given compression term specified in the "
                         "MaterialModel output PrescribedPlasticDilation.");
    }
    prm.leave_subsection();

    // next declare parameters that pertain to the equations to be
    // solved, along with boundary conditions etc. note that at this
    // point we do not know yet which geometry model we will use, so
    // we do not know which symbolic names will be valid to address individual
    // parts of the boundary. we can only work around this by allowing any string
    // to indicate a boundary
    prm.enter_subsection ("Melt settings");
    {
      prm.declare_entry ("Include melt transport", "false",
                         Patterns::Bool (),
                         "Whether to include the transport of melt into the model or not. If this "
                         "is set to true, two additional pressures (the fluid pressure and the "
                         "compaction pressure) will be added to the finite element. "
                         "Including melt transport in the simulation also requires that there is "
                         "one compositional field that has the name `porosity'. This field will "
                         "be used for computing the additional pressures and the melt velocity, "
                         "and has a different advection equation than other compositional fields, "
                         "as it is effectively advected with the melt velocity.");
    }
    prm.leave_subsection();

    prm.enter_subsection ("Mesh deformation");
    {
      prm.declare_entry ("Mesh deformation boundary indicators", "",
                         Patterns::List (Patterns::Anything()),
                         "A comma separated list of names denoting those boundaries "
                         "where there is some type of mesh deformation. Set to nothing to disable all "
                         "deformation computations."
                         "\n\n"
                         "The names of the boundaries listed here can either by "
                         "numbers (in which case they correspond to the numerical "
                         "boundary indicators assigned by the geometry object), or they "
                         "can correspond to any of the symbolic names the geometry object "
                         "may have provided for each part of the boundary. You may want "
                         "to compare this with the documentation of the geometry model you "
                         "use in your model.");
    }
    prm.leave_subsection();

    prm.enter_subsection ("Boundary traction model");
    {
      prm.declare_entry ("Prescribed traction boundary indicators", "",
                         Patterns::Map (Patterns::Anything(),
                                        Patterns::Selection(BoundaryTraction::get_names<dim>())),
                         "A comma separated list denoting those boundaries "
                         "on which a traction force is prescribed, i.e., where "
                         "known external forces act, resulting in an unknown velocity. This is "
                         "often used to model ``open'' boundaries where we only know the pressure. "
                         "This pressure then produces a force that is normal to the boundary and "
                         "proportional to the pressure."
                         "\n\n"
                         "The format of valid entries for this parameter is that of a map "
                         "given as ``key1 [selector]: value1, key2 [selector]: value2, key3: value3, ...'' where "
                         "each key must be a valid boundary indicator (which is either an "
                         "integer or the symbolic name the geometry model in use may have "
                         "provided for this part of the boundary) "
                         "and each value must be one of the currently implemented boundary "
                         "traction models. ``selector'' is an optional string given as a subset "
                         "of the letters `xyz' that allows you to apply the boundary conditions "
                         "only to the components listed. As an example, '1 y: function' applies "
                         "the type `function' to the y component on boundary 1. Without a selector "
                         "it will affect all components of the traction.");
    }
    prm.leave_subsection();

    prm.enter_subsection ("Boundary heat flux model");
    {
      prm.declare_entry ("Fixed heat flux boundary indicators", "",
                         Patterns::List (Patterns::Anything()),
                         "A comma separated list of names denoting those boundaries "
                         "on which the heat flux is fixed and described by the "
                         "boundary heat flux object selected in the 'Model name' parameter. "
                         "All boundary indicators used by the geometry but not explicitly "
                         "listed here or in the list of 'Fixed temperature boundary indicators' "
                         "in the 'Boundary temperature model' will end up with no-flux "
                         "(insulating) boundary conditions."
                         "\n\n"
                         "The names of the boundaries listed here can either be "
                         "numbers (in which case they correspond to the numerical "
                         "boundary indicators assigned by the geometry object), or they "
                         "can correspond to any of the symbolic names the geometry object "
                         "may have provided for each part of the boundary. You may want "
                         "to compare this with the documentation of the geometry model you "
                         "use in your model."
                         "\n\n"
                         "This parameter only describes which boundaries have a fixed "
                         "heat flux, but not what heat flux should hold on these "
                         "boundaries. The latter piece of information needs to be "
                         "implemented in a plugin in the BoundaryHeatFlux "
                         "group, unless an existing implementation in this group "
                         "already provides what you want.");
    }
    prm.leave_subsection();


    prm.enter_subsection ("Nullspace removal");
    {
      prm.declare_entry ("Remove nullspace", "",
                         Patterns::MultipleSelection("net rotation|angular momentum|"
                                                     "net translation|linear momentum|"
                                                     "net x translation|net y translation|net z translation|"
                                                     "linear x momentum|linear y momentum|linear z momentum"),
                         "Choose none, one or several from "
                         "\n\n"
                         "\\begin{itemize} \\item net rotation \\item angular momentum \\item net translation "
                         "\\item linear momentum \\item net x translation \\item net y translation "
                         "\\item net z translation \\item linear x momentum \\item linear y momentum "
                         "\\item linear z momentum \\end{itemize}"
                         "\n\n"
                         "These are a selection of operations to remove certain parts of the nullspace from "
                         "the velocity after solving. For some geometries and certain boundary conditions "
                         "the velocity field is not uniquely determined but contains free translations "
                         "and/or rotations. Depending on what you specify here, these non-determined "
                         "modes will be removed from the velocity field at the end of the Stokes solve step.\n"
                         "\n\n"
                         "The ``angular momentum'' option removes a rotation such that the net angular momentum "
                         "is zero. The ``linear * momentum'' options remove translations such that the net "
                         "momentum in the relevant direction is zero.  The ``net rotation'' option removes the "
                         "net rotation of the domain, and the ``net * translation'' options remove the "
                         "net translations in the relevant directions.  For most problems there should not be a "
                         "significant difference between the momentum and rotation/translation versions of "
                         "nullspace removal, although the momentum versions are more physically motivated. "
                         "They are equivalent for constant density simulations, and approximately equivalent "
                         "when the density variations are small."
                         "\n\n"
                         "Note that while more than one operation can be selected it only makes sense to "
                         "pick one rotational and one translational operation.");
    }
    prm.leave_subsection();

    prm.enter_subsection ("Mesh refinement");
    {
      prm.declare_entry ("Initial global refinement", "2",
                         Patterns::Integer (0),
                         "The number of global refinement steps performed on "
                         "the initial coarse mesh, before the problem is first "
                         "solved there.");
      prm.declare_entry ("Initial adaptive refinement", "0",
                         Patterns::Integer (0),
                         "The number of adaptive refinement steps performed after "
                         "initial global refinement but while still within the first "
                         "time step. These refinement steps (n) are added to the value "
                         "for initial global refinement (m) so that the final mesh has "
                         "cells that are at most on refinement level $n+m$.");
      prm.declare_entry ("Time steps between mesh refinement", "10",
                         Patterns::Integer (0),
                         "The number of time steps after which the mesh is to be "
                         "adapted again based on computed error indicators. If 0 "
                         "then the mesh will never be changed.");
      prm.declare_entry ("Refinement fraction", "0.3",
                         Patterns::Double(0., 1.),
                         "Cells are sorted from largest to smallest by their total error "
                         "(determined by the Strategy). Then the cells with the largest "
                         "error (top of this sorted list) that account for given fraction "
                         "of the error are refined.");
      prm.declare_entry ("Coarsening fraction", "0.05",
                         Patterns::Double(0., 1.),
                         "Cells are sorted from largest to smallest by their total error "
                         "(determined by the Strategy). Then the cells with the smallest "
                         "error (bottom of this sorted list) that account for the given fraction "
                         "of the error are coarsened.");
      prm.declare_entry ("Adapt by fraction of cells", "false",
                         Patterns::Bool(),
                         "Use fraction of the total number of cells instead of "
                         "fraction of the total error as the limit for refinement "
                         "and coarsening.");
      prm.declare_entry ("Minimum refinement level", "0",
                         Patterns::Integer (0),
                         "The minimum refinement level each cell should have, "
                         "and that can not be exceeded by coarsening. "
                         "Should not be higher than the 'Initial global refinement' "
                         "parameter.");
      prm.declare_entry ("Additional refinement times", "",
                         Patterns::List (Patterns::Double (0.)),
                         "A list of times so that if the end time of a time step "
                         "is beyond this time, an additional round of mesh refinement "
                         "is triggered. This is mostly useful to make sure we "
                         "can get through the initial transient phase of a simulation "
                         "on a relatively coarse mesh, and then refine again when we "
                         "are in a time range that we are interested in and where "
                         "we would like to use a finer mesh. Units: Each element of the "
                         "list has units years if the "
                         "'Use years in output instead of seconds' parameter is set; "
                         "seconds otherwise.");
      prm.declare_entry ("Run postprocessors on initial refinement", "false",
                         Patterns::Bool (),
                         "Whether or not the postprocessors should be executed after "
                         "each of the initial adaptive refinement cycles that are run at "
                         "the start of the simulation. This is useful for "
                         "plotting/analyzing how the mesh refinement parameters are "
                         "working for a particular model.");
      prm.declare_entry ("Skip solvers on initial refinement", "false",
                         Patterns::Bool (),
                         "Whether or not solvers should be executed during the initial "
                         "adaptive refinement cycles that are run at the start of the "
                         "simulation.");
      prm.declare_entry ("Skip setup initial conditions on initial refinement", "false",
                         Patterns::Bool (),
                         "Whether or not the initial conditions should be set up during the "
                         "adaptive refinement cycles that are run at the start of the "
                         "simulation.");
    }
    prm.leave_subsection();

    prm.enter_subsection ("Postprocess");
    {
      prm.declare_entry ("Run postprocessors on nonlinear iterations", "false",
                         Patterns::Bool (),
                         "Whether or not the postprocessors should be executed after "
                         "each of the nonlinear iterations done within one time step. "
                         "As this is mainly an option for the purposes of debugging, "
                         "it is not supported when the 'Time between graphical output' "
                         "is larger than zero, or when the postprocessor is not intended "
                         "to be run more than once per timestep.");
    }
    prm.leave_subsection();

    prm.enter_subsection ("Checkpointing");
    {
      prm.declare_entry ("Time between checkpoint", "0",
                         Patterns::Integer (0),
                         "The wall time between performing checkpoints. "
                         "If 0, will use the checkpoint step frequency instead. "
                         "Units: Seconds.");
      prm.declare_entry ("Steps between checkpoint", "0",
                         Patterns::Integer (0),
                         "The number of timesteps between performing checkpoints. "
                         "If 0 and time between checkpoint is not specified, "
                         "checkpointing will not be performed. "
                         "Units: None.");
    }
    prm.leave_subsection ();

    prm.enter_subsection ("Discretization");
    {
      prm.declare_entry ("Stokes velocity polynomial degree", "2",
                         Patterns::Integer (1),
                         "The polynomial degree to use for the velocity variables "
                         "in the Stokes system. The polynomial degree for the pressure "
                         "variable will then be one less in order to make the velocity/pressure "
                         "pair conform with the usual LBB (Babu{\\v s}ka-Brezzi) condition. In "
                         "other words, we are using a Taylor-Hood element for the Stokes "
                         "equations and this parameter indicates the polynomial degree of it. "
                         "As an example, a value of 2 for this parameter will yield the "
                         "element $Q_2^d \\times Q_1$ for the $d$ velocity components and the "
                         "pressure, respectively (unless the `Use locally conservative "
                         "discretization' parameter is set, which modifies the pressure "
                         "element). "
                         "\n\n"
                         "Be careful if you choose 1 as the degree. The resulting element "
                         "is not stable and it may lead to artifacts in the solution. "
                         "Units: None.");
      prm.declare_entry ("Temperature polynomial degree", "2",
                         Patterns::Integer (1),
                         "The polynomial degree to use for the temperature variable. "
                         "As an example, a value of 2 for this parameter will yield "
                         "either the element $Q_2$ or $DGQ_2$ for the temperature "
                         "field, depending on whether we use a continuous or "
                         "discontinuous field. "
                         "Units: None.");
      prm.declare_entry ("Composition polynomial degree", "2",
                         Patterns::Integer (0),
                         "The polynomial degree to use for the composition variable(s). "
                         "As an example, a value of 2 for this parameter will yield "
                         "either the element $Q_2$ or $DGQ_2$ for the compositional "
                         "field(s), depending on whether we use continuous or "
                         "discontinuous field(s). "
                         "\n\n"
                         "For continuous elements, the value needs to be 1 or larger "
                         "as $Q_1$ is the lowest order element, while $DGQ_0$ is a "
                         "valid choice. "
                         "Units: None.");
      prm.declare_entry ("Use locally conservative discretization", "false",
                         Patterns::Bool (),
                         "Whether to use a Stokes discretization that is locally "
                         "conservative at the expense of a larger number of degrees "
                         "of freedom (true), or to go with a cheaper discretization "
                         "that does not locally conserve mass, although it is "
                         "globally conservative (false)."
                         "\n\n"
                         "When using a locally "
                         "conservative discretization, the finite element space for "
                         "the pressure is discontinuous between cells and is the "
                         "polynomial space $P_{-(k-1)}$ of polynomials of degree $k-1$ in "
                         "each variable separately. Here, $k$ is the value "
                         "given in the parameter ``Stokes velocity polynomial degree'', and "
                         "consequently the polynomial degree for the pressure, $k-1$, is one "
                         "lower than that for the velocity."
                         "\n\n"
                         "As a consequence of choosing this element for the pressure rather "
                         "than the more commonly used $Q_{k-1}$ element that is continuous, "
                         "it can be shown that if the medium is considered incompressible "
                         "then the computed discrete velocity "
                         "field $\\mathbf u_h$ satisfies the property $\\int_ {\\partial K} \\mathbf u_h "
                         "\\cdot \\mathbf n = 0$ for every cell $K$, i.e., for each cell inflow and "
                         "outflow exactly balance each other as one would expect for an "
                         "incompressible medium. In other words, the velocity field is "
                         "\\textit{locally conservative}."
                         "\n\n"
                         "On the other hand, if this parameter is set to ``false''"
                         "(the default), then the finite element space is chosen as $Q_{k-1}$. "
                         "This choice does not yield the local conservation property but "
                         "has the advantage of requiring fewer degrees of freedom. Furthermore, "
                         "the error is generally smaller with this choice."
                         "\n\n"
                         "For an in-depth discussion of these issues and a quantitative evaluation "
                         "of the different choices, see \\cite{KHB12}.");
      prm.declare_entry ("Use equal order interpolation for Stokes", "false",
                         Patterns::Bool(),
                         "By default (i.e., when this parameter is set to its default value "
                         "`false') \\aspect{} uses finite element combinations in which the "
                         "pressure shape functions are polynomials one degree lower than "
                         "the shape functions for the velocity. An example is the "
                         "Taylor-Hood element that uses $Q_k$ elements for the velocity "
                         "and $Q_{k-1}$ for the pressure. This is because using the "
                         "\\textit{same} polynomial degree for both the velocity and the "
                         "pressure turns out to violate some mathematical properties "
                         "necessary to make the problem solvable. (In particular, the"
                         "condition in question goes by the name ``inf-sup'' or "
                         "Babu{\\v s}ka-Brezzi or LBB condition.) A consequence of "
                         "violating this condition is that the pressure may show "
                         "oscillations and not converge to the correct pressure."
                         "\n\n"
                         "That said, people have often used $Q_1$ elements for both the"
                         "velocity and pressure anyway. This is commonly referred to as "
                         "using the Q1-Q1 method. It is, by default, not stable as "
                         "mentioned above, but it can be made stable by adding small "
                         "amount of compressibility to the model. There are numerous "
                         "ways to do that. Today, the way that is generally considered "
                         "to be the best approach is the one by Dohrmann and Bochev "
                         "\\cite{DohrmannBochev2004}."
                         "\n\n"
                         "When this parameter is set to ``true'', then \\aspect{} "
                         "will use this method by using $Q_k\times Q_k$ elements for "
                         "velocity and pressure, respectively, where $k$ is the value "
                         "provided for the parameter ``Stokes velocity polynomial "
                         "degree''."
                         "\n\n"
                         "\\note{While \\aspect{} \\textit{allows} you to use this "
                         "  method, it is generally understood that this is not a "
                         "  great idea as it leads to rather low accuracy in "
                         "  general. It also leads to substantial problems when "
                         "  using free surfaces. As a consequence, the presence "
                         "  of this parameter should not be seen as an "
                         "  endorsement of the method, or a suggestion to "
                         "  actually use it. It simply makes the method available.}");
      prm.declare_entry ("Use discontinuous temperature discretization", "false",
                         Patterns::Bool (),
                         "Whether to use a temperature discretization that is discontinuous "
                         "as opposed to continuous. This then requires the assembly of face terms "
                         "between cells, and weak imposition of boundary terms for the temperature "
                         "field via the interior-penalty discontinuous Galerkin method.");
      prm.declare_entry ("Use discontinuous composition discretization", "false",
                         Patterns::Bool (),
                         "Whether to use a composition discretization that is discontinuous "
                         "as opposed to continuous. This then requires the assembly of face terms "
                         "between cells, and weak imposition of boundary terms for the composition "
                         "field via the discontinuous Galerkin method.");

      prm.enter_subsection ("Stabilization parameters");
      {
        prm.declare_entry ("Stabilization method", "entropy viscosity",
                           Patterns::Selection("entropy viscosity|SUPG"),
                           "Select the method for stabilizing the advection equation. The original "
                           "method implemented is 'entropy viscosity' as described in \\cite {KHB12}. "
                           "SUPG is currently experimental.");

        prm.declare_entry ("Use artificial viscosity smoothing", "false",
                           Patterns::Bool (),
                           "If set to false, the artificial viscosity of a cell is computed and "
                           "is computed on every cell separately as discussed in \\cite{KHB12}. "
                           "If set to true, the maximum of the artificial viscosity in "
                           "the cell as well as the neighbors of the cell is computed and used "
                           "instead.");

        prm.declare_entry ("alpha", "2",
                           Patterns::Integer (1, 2),
                           "The exponent $\\alpha$ in the entropy viscosity stabilization. Valid "
                           "options are 1 or 2. The recommended setting is 2. (This parameter does "
                           "not correspond to any variable in the 2012 paper by Kronbichler, "
                           "Heister and Bangerth that describes ASPECT, see \\cite{KHB12}. "
                           "Rather, the paper always uses 2 as the exponent in the definition "
                           "of the entropy, following equation (15) of the paper. The full "
                           "approach is discussed in \\cite{GPP11}.) Note that this is not the "
                           "thermal expansion coefficient, also commonly referred to as $\\alpha$."
                           "Units: None.");
        prm.declare_entry ("cR", "0.11",
                           Patterns::List(Patterns::Double (0.)),
                           "The $c_R$ factor in the entropy viscosity "
                           "stabilization. This parameter controls the part of the entropy viscosity "
                           "that depends on the solution field itself and its residual in addition "
                           "to the cell diameter and the maximum velocity in the cell. "
                           "This parameter can be given as a single value or as a list with as "
                           "many entries as one plus the number of compositional fields. In the "
                           "former case all advection fields use the same stabilization parameters, "
                           "in the latter case each field (temperature first, then all compositions) "
                           "use individual parameters. This can be useful to reduce the stabilization "
                           "for the temperature, which already has some physical diffusion. "
                           "(For historical reasons, the name used here is different "
                           "from the one used in the 2012 paper by Kronbichler, "
                           "Heister and Bangerth that describes ASPECT, see \\cite{KHB12}. "
                           "This parameter corresponds "
                           "to the factor $\\alpha_E$ in the formulas following equation (15) of "
                           "the paper.) Units: None.");
        prm.declare_entry ("beta", "0.052",
                           Patterns::List(Patterns::Double (0.)),
                           "The $\\beta$ factor in the artificial viscosity "
                           "stabilization. This parameter controls the maximum dissipation of the "
                           "entropy viscosity, which is the part that only scales with the cell diameter "
                           "and the maximum velocity in the cell, but does not depend on the solution "
                           "field itself or its residual. An appropriate value for 2d is 0.052 and "
                           "0.78 for 3d. (For historical reasons, the name used here is different "
                           "from the one used in the 2012 paper by Kronbichler, "
                           "Heister and Bangerth that describes ASPECT, see \\cite{KHB12}. "
                           "This parameter can be given as a single value or as a list with as "
                           "many entries as one plus the number of compositional fields. In the "
                           "former case all advection fields use the same stabilization parameters, "
                           "in the latter case each field (temperature first, then all compositions) "
                           "use individual parameters. This can be useful to reduce the stabilization "
                           "for the temperature, which already has some physical diffusion. "
                           "This parameter corresponds "
                           "to the factor $\\alpha_{\\text{max}}$ in the formulas following equation (15) of "
                           "the paper.) Units: None.");
        prm.declare_entry ("gamma", "0.0",
                           Patterns::Double (0.),
                           "The strain rate scaling factor in the artificial viscosity "
                           "stabilization. This parameter determines how much the strain rate (in addition "
                           "to the velocity) should influence the stabilization. (This parameter does "
                           "not correspond to any variable in the 2012 paper by Kronbichler, "
                           "Heister and Bangerth that describes ASPECT, see \\cite{KHB12}. "
                           "Rather, the paper always uses "
                           "0, i.e. they specify the maximum dissipation $\\nu_h^\\text{max}$ as "
                           "$\\nu_h^\\text{max}\\vert_K = \\alpha_{\\text{max}} h_K \\|\\mathbf u\\|_{\\infty,K}$. "
                           "Here, we use "
                           "$\\|\\lvert\\mathbf u\\rvert + \\gamma h_K \\lvert\\varepsilon (\\mathbf u)\\rvert\\|_{\\infty,K}$ "
                           "instead of $\\|\\mathbf u\\|_{\\infty,K}$. "
                           "Units: None.");
        prm.declare_entry ("Discontinuous penalty", "10.",
                           Patterns::Double (0.),
                           "The value used to penalize discontinuities in the discontinuous Galerkin "
                           "method. This is used only for the temperature field, and not for the composition "
                           "field, as pure advection does not use the interior penalty method. This "
                           "is largely empirically decided -- it must be large enough to ensure "
                           "the bilinear form is coercive, but not so large as to penalize "
                           "discontinuity at all costs.");
        prm.declare_entry ("Use limiter for discontinuous temperature solution", "false",
                           Patterns::Bool (),
                           "Whether to apply the bound preserving limiter as a correction after computing "
                           "the discontinuous temperature solution. Currently we apply this only to the "
                           "temperature solution if the 'Global temperature maximum' and "
                           "'Global temperature minimum' are already defined in the .prm file. "
                           "This limiter keeps the discontinuous solution in the range given by "
                           "'Global temperature maximum' and 'Global temperature minimum'.");
        prm.declare_entry ("Use limiter for discontinuous composition solution", "false",
                           Patterns::Bool (),
                           "Whether to apply the bound preserving limiter as a correction after having "
                           "the discontinuous composition solution. Currently we apply this only to the "
                           "compositional solution if the 'Global composition maximum' and "
                           "'Global composition minimum' are already defined in the .prm file. "
                           "This limiter keeps the discontinuous solution in the range given by "
                           "Global composition maximum' and 'Global composition minimum'.");
        prm.declare_entry ("Global temperature maximum",
                           boost::lexical_cast<std::string>(std::numeric_limits<double>::max()),
                           Patterns::Double (),
                           "The maximum global temperature value that will be used in the bound preserving "
                           "limiter for the discontinuous solutions from temperature advection fields.");
        prm.declare_entry ("Global temperature minimum",
                           boost::lexical_cast<std::string>(-std::numeric_limits<double>::max()),
                           Patterns::Double (),
                           "The minimum global temperature value that will be used in the bound preserving "
                           "limiter for the discontinuous solutions from temperature advection fields.");
        prm.declare_entry ("Global composition maximum",
                           boost::lexical_cast<std::string>(std::numeric_limits<double>::max()),
                           Patterns::List(Patterns::Double ()),
                           "The maximum global composition values that will be used in the bound preserving "
                           "limiter for the discontinuous solutions from composition advection fields. "
                           "The number of the input 'Global composition maximum' values separated by ',' has to be "
                           "one or the same as the number of the compositional fields. When only one value "
                           "is supplied, this same value is assumed for all compositional fields.");
        prm.declare_entry ("Global composition minimum",
                           boost::lexical_cast<std::string>(-std::numeric_limits<double>::max()),
                           Patterns::List(Patterns::Double ()),
                           "The minimum global composition value that will be used in the bound preserving "
                           "limiter for the discontinuous solutions from composition advection fields. "
                           "The number of the input 'Global composition minimum' values separated by ',' has to be "
                           "one or the same as the number of the compositional fields. When only one value "
                           "is supplied, this same value is assumed for all compositional fields.");
      }
      prm.leave_subsection ();
    }
    prm.leave_subsection ();

    prm.enter_subsection ("Temperature field");
    {
      prm.declare_entry ("Temperature method", "field",
                         Patterns::Selection("field|prescribed field"),
                         "A comma separated list denoting the solution method of the "
                         "temperature field. Each entry of the list must be "
                         "one of the currently implemented field types."
                         "\n\n"
                         "These choices correspond to the following methods by which "
                         "the temperature field gains its values:"
                         "\\begin{itemize}"
                         "\\item ``field'': If the temperature is marked with this "
                         "method, then its values are computed in each time step by "
                         "solving the temperature advection-diffusion equation. In other words, "
                         "this corresponds to the usual notion of a temperature. "
                         "\n"
                         "\\item ``prescribed field'': The value of the temperature is determined "
                         "in each time step from the material model. If a compositional field is "
                         "marked with this method, then the value of a specific additional material "
                         "model output, called the `PrescribedTemperatureOutputs' is interpolated "
                         "onto the temperature. This field does not change otherwise, it is not "
                         "advected with the flow. \n"
                         "\\end{itemize}");
    }
    prm.leave_subsection();

    prm.enter_subsection ("Compositional fields");
    {
      prm.declare_entry ("Number of fields", "0",
                         Patterns::Integer (0),
                         "The number of fields that will be advected along with the flow field, excluding "
                         "velocity, pressure and temperature.");
      prm.declare_entry ("Names of fields", "",
                         Patterns::List(Patterns::Anything()),
                         "A user-defined name for each of the compositional fields requested.");
      prm.declare_entry ("Compositional field methods", "",
                         Patterns::List (Patterns::Selection("field|particles|volume of fluid|static|melt field|prescribed field|prescribed field with diffusion")),
                         "A comma separated list denoting the solution method of each "
                         "compositional field. Each entry of the list must be "
                         "one of the currently implemented field types."
                         "\n\n"
                         "These choices correspond to the following methods by which "
                         "compositional fields gain their values:"
                         "\\begin{itemize}"
                         "\\item ``field'': If a compositional field is marked with this "
                         "method, then its values are computed in each time step by "
                         "advecting along the values of the previous time step using the "
                         "velocity field, and applying reaction rates to it. In other words, "
                         "this corresponds to the usual notion of a composition field as "
                         "mentioned in Section~\\ref{sec:compositional}. "
                         "\n"
                         "\\item ``particles'': If a compositional field is marked with "
                         "this method, then its values are obtained in each time step "
                         "by interpolating the corresponding properties from the "
                         "particles located on each cell. The time evolution therefore "
                         "happens because particles move along with the velocity field, "
                         "and particle properties can react with each other as well. "
                         "See Section~\\ref{sec:particles} for more information about "
                         "how particles behave."
                         "\n"
                         "\\item ``volume of fluid``: If a compositional field "
                         "is marked with this method, then its values are "
                         "obtained in each timestep by reconstructing a "
                         "polynomial finite element approximation on each cell "
                         "from a volume of fluid interface tracking method, "
                         "which is used to compute the advection updates."
                         "\n"
                         "\\item ``static'': If a compositional field is marked "
                         "this way, then it does not evolve at all. Its values are "
                         "simply set to the initial conditions, and will then "
                         "never change."
                         "\n"
                         "\\item ``melt field'': If a compositional field is marked with this "
                         "method, then its values are computed in each time step by "
                         "advecting along the values of the previous time step using the "
                         "melt velocity, and applying reaction rates to it. In other words, "
                         "this corresponds to the usual notion of a composition field as "
                         "mentioned in Section~\\ref{sec:compositional}, except that it is "
                         "advected with the melt velocity instead of the solid velocity. "
                         "This method can only be chosen if melt transport is active in the "
                         "model."
                         "\n"
                         "\\item ``prescribed field'': The value of these fields is determined "
                         "in each time step from the material model. If a compositional field is "
                         "marked with this method, then the value of a specific additional material "
                         "model output, called the `PrescribedFieldOutputs' is interpolated "
                         "onto the field. This field does not change otherwise, it is not "
                         "advected with the flow."
                         "\n"
                         "\\item ``prescribed field with diffusion'': If a compositional field is "
                         "marked this way, the value of a specific additional material model output, "
                         "called the `PrescribedFieldOutputs' is interpolated onto the field, as in "
                         "the ``prescribed field'' method. Afterwards, the field is diffused based on "
                         "a solver parameter, the diffusion length scale, smoothing the field. "
                         "Specifically, the field is updated by solving the equation "
                         "$(I-l^2 \\Delta) C_\\text{smoothed} = C_\\text{prescribed}$, "
                         "where $l$ is the diffusion length scale. Note that this means that the amount "
                         "of diffusion is independent of the time step size, and that the field is not "
                         "advected with the flow."
                         "\\end{itemize}");
      prm.declare_entry ("Mapped particle properties", "",
                         Patterns::Map (Patterns::Anything(),
                                        Patterns::Anything()),
                         "A comma separated list denoting the particle properties "
                         "that will be projected to those compositional fields that "
                         "are of the ``particles'' field type."
                         "\n\n"
                         "The format of valid entries for this parameter is that of a map "
                         "given as ``key1: value1, key2: value2 [component2], key3: value3 [component4], "
                         "...'' where each key must be a valid field name of the "
                         "``particles'' type, and each value must be one of the currently "
                         "selected particle properties. Component is a component index of "
                         "the particle property that is 0 by default, but can be set up to "
                         "n-1, where n is the number of vector components of this particle "
                         "property. The component indicator only needs to be "
                         "set if not the first component of the particle property should be "
                         "mapped (e.g. the $y$-component of the velocity at the particle positions).");
      prm.declare_entry ("List of normalized fields", "",
                         Patterns::List (Patterns::Integer(0)),
                         "A list of integers smaller than or equal to the number of "
                         "compositional fields. All compositional fields in this "
                         "list will be normalized before the first timestep. "
                         "The normalization is implemented in the following way: "
                         "First, the sum of the fields to be normalized is calculated "
                         "at every point and the global maximum is determined. "
                         "Second, the compositional fields to be normalized are "
                         "divided by this maximum.");
    }
    prm.leave_subsection ();

    // Finally declare a couple of parameters related how we should
    // evaluate the material models when assembling the matrix and
    // preconditioner
    prm.enter_subsection ("Material model");
    {
      prm.declare_entry ("Material averaging", "none",
                         Patterns::Selection(MaterialModel::MaterialAveraging::
                                             get_averaging_operation_names()),
                         "Whether or not (and in the first case, how) to do any averaging of "
                         "material model output data when constructing the linear systems "
                         "for velocity/pressure, temperature, and compositions in each "
                         "time step, as well as their corresponding preconditioners."
                         "\n\n"
                         "Possible choices: " + MaterialModel::MaterialAveraging::
                         get_averaging_operation_names()
                         +
                         "\n\n"
                         "The process of averaging, and where it may be used, is "
                         "discussed in more detail in "
                         "Section~\\ref{sec:sinker-with-averaging}."
                         "\n\n"
                         "More averaging schemes are available in the averaging material "
                         "model. This material model is a ``compositing material model'' "
                         "which can be used in combination with other material models.");
    }
    prm.leave_subsection ();

    prm.enter_subsection ("Volume of Fluid");
    {
      prm.declare_entry ("Enable interface tracking", "false",
                         Patterns::Bool (),
                         "When set to true, Volume of Fluid interface tracking will be used");
    }
    prm.leave_subsection ();

    // declare the VolumeOfFluid parameters
    VolumeOfFluidHandler<dim>::declare_parameters(prm);

    // then, finally, let user additions that do not go through the usual
    // plugin mechanism, declare their parameters if they have subscribed
    // to the relevant signals
    SimulatorSignals<dim>::declare_additional_parameters (dim, prm);
  }



  template <int dim>
  void
  Parameters<dim>::
  parse_parameters (ParameterHandler &prm,
                    const MPI_Comm mpi_communicator)
  {
    // first, make sure that the ParameterHandler parser agrees
    // with the code in main() about the meaning of the "Dimension"
    // parameter
    AssertThrow (prm.get_integer("Dimension") == dim,
                 ExcInternalError());

    CFL_number              = prm.get_double ("CFL number");
    use_conduction_timestep = prm.get_bool ("Use conduction timestep");
    convert_to_years        = prm.get_bool ("Use years in output instead of seconds");
    timing_output_frequency = prm.get_integer ("Timing output frequency");
    world_builder_file      = prm.get("World builder file");

    maximum_time_step       = prm.get_double("Maximum time step");
    if (convert_to_years == true)
      maximum_time_step *= year_in_seconds;

    maximum_relative_increase_time_step = prm.get_double("Maximum relative increase in time step") * 0.01;
    maximum_first_time_step = prm.get_double("Maximum first time step");
    if (convert_to_years == true)
      maximum_first_time_step *= year_in_seconds;

    {
      const std::string solver_scheme = prm.get ("Nonlinear solver scheme");
      if (solver_scheme == "single Advection, single Stokes" || solver_scheme == "IMPES")
        nonlinear_solver = NonlinearSolver::single_Advection_single_Stokes;
      else if (solver_scheme == "iterated Advection and Stokes" || solver_scheme == "iterated IMPES")
        nonlinear_solver = NonlinearSolver::iterated_Advection_and_Stokes;
      else if (solver_scheme == "single Advection, iterated Stokes" || solver_scheme == "iterated Stokes")
        nonlinear_solver = NonlinearSolver::single_Advection_iterated_Stokes;
      else if (solver_scheme == "no Advection, iterated Stokes" || solver_scheme == "Stokes only")
        nonlinear_solver = NonlinearSolver::no_Advection_iterated_Stokes;
      else if (solver_scheme == "no Advection, single Stokes")
        nonlinear_solver = NonlinearSolver::no_Advection_single_Stokes;
      else if (solver_scheme == "no Advection, iterated defect correction Stokes")
        nonlinear_solver = NonlinearSolver::no_Advection_iterated_defect_correction_Stokes;
      else if (solver_scheme == "single Advection, iterated defect correction Stokes")
        nonlinear_solver = NonlinearSolver::single_Advection_iterated_defect_correction_Stokes;
      else if (solver_scheme == "iterated Advection and defect correction Stokes")
        nonlinear_solver = NonlinearSolver::iterated_Advection_and_defect_correction_Stokes;
      else if (solver_scheme == "iterated Advection and Newton Stokes" || solver_scheme == "Newton Stokes")
        nonlinear_solver = NonlinearSolver::iterated_Advection_and_Newton_Stokes;
      else if (solver_scheme == "single Advection, iterated Newton Stokes")
        nonlinear_solver = NonlinearSolver::single_Advection_iterated_Newton_Stokes;
      else if (solver_scheme == "single Advection, no Stokes" || solver_scheme == "Advection only")
        nonlinear_solver = NonlinearSolver::single_Advection_no_Stokes;
      else if (solver_scheme == "first timestep only, single Stokes")
        nonlinear_solver = NonlinearSolver::first_timestep_only_single_Stokes;
      else if (solver_scheme == "no Advection, no Stokes")
        nonlinear_solver = NonlinearSolver::no_Advection_no_Stokes;
      else
        AssertThrow (false, ExcNotImplemented());
    }

    prm.enter_subsection ("Solver parameters");
    {
      temperature_solver_tolerance    = prm.get_double ("Temperature solver tolerance");
      composition_solver_tolerance    = prm.get_double ("Composition solver tolerance");

      prm.enter_subsection ("Advection solver parameters");
      {
        advection_gmres_restart_length     = prm.get_integer("GMRES solver restart length");
      }
      prm.leave_subsection ();

      prm.enter_subsection ("Stokes solver parameters");
      {
        stokes_solver_type = StokesSolverType::parse(prm.get("Stokes solver type"));
        if (prm.get_bool("Use direct solver for Stokes system"))
          stokes_solver_type = StokesSolverType::direct_solver;
        use_direct_stokes_solver        = stokes_solver_type==StokesSolverType::direct_solver;
        stokes_krylov_type = StokesKrylovType::parse(prm.get("Krylov method for cheap solver steps"));
        idr_s_parameter    = prm.get_integer("IDR(s) parameter");

        linear_stokes_solver_tolerance  = prm.get_double ("Linear solver tolerance");
        n_cheap_stokes_solver_steps     = prm.get_integer ("Number of cheap Stokes solver steps");
        n_expensive_stokes_solver_steps = prm.get_integer ("Maximum number of expensive Stokes solver steps");
        linear_solver_A_block_tolerance = prm.get_double ("Linear solver A block tolerance");
        use_full_A_block_preconditioner = prm.get_bool ("Use full A block as preconditioner");
        linear_solver_S_block_tolerance = prm.get_double ("Linear solver S block tolerance");
        stokes_gmres_restart_length     = prm.get_integer("GMRES solver restart length");
      }
      prm.leave_subsection ();

      prm.enter_subsection ("AMG parameters");
      {
        AMG_smoother_type                      = prm.get ("AMG smoother type");
        AMG_smoother_sweeps                    = prm.get_integer ("AMG smoother sweeps");
        AMG_aggregation_threshold              = prm.get_double ("AMG aggregation threshold");
        AMG_output_details                     = prm.get_bool ("AMG output details");
      }
      prm.leave_subsection ();
      prm.enter_subsection ("Operator splitting parameters");
      {
        reaction_time_step       = prm.get_double("Reaction time step");
        AssertThrow (reaction_time_step > 0,
                     ExcMessage("Reaction time step must be greater than 0."));
        if (convert_to_years == true)
          reaction_time_step *= year_in_seconds;
        reaction_steps_per_advection_step = prm.get_integer ("Reaction time steps per advection step");
      }
      prm.leave_subsection ();
      prm.enter_subsection ("Diffusion solver parameters");
      {
        diffusion_length_scale = prm.get_double("Diffusion length scale");
      }
      prm.leave_subsection ();
    }
    prm.leave_subsection ();

    nonlinear_tolerance = prm.get_double("Nonlinear solver tolerance");

    max_nonlinear_iterations = prm.get_integer ("Max nonlinear iterations");
    max_nonlinear_iterations_in_prerefinement = prm.get_integer ("Max nonlinear iterations in pre-refinement");

    start_time              = prm.get_double ("Start time");
    if (convert_to_years == true)
      start_time *= year_in_seconds;

    output_directory        = prm.get ("Output directory");
    if (output_directory.size() == 0)
      output_directory = "./";
    else if (output_directory[output_directory.size()-1] != '/')
      output_directory += "/";

    Utilities::create_directory (output_directory,
                                 mpi_communicator,
                                 false);

    if (prm.get ("Resume computation") == "true")
      resume_computation = true;
    else if (prm.get ("Resume computation") == "false")
      resume_computation = false;
    else if (prm.get ("Resume computation") == "auto")
      {
        resume_computation = Utilities::fexists(output_directory+"restart.mesh");
      }
    else
      AssertThrow (false, ExcMessage ("Resume computation parameter must be either `true', `false', or `auto'."));
#ifndef DEAL_II_WITH_ZLIB
    AssertThrow (resume_computation == false,
                 ExcMessage ("You need to have deal.II configured with the `libz' "
                             "option if you want to resume a computation from a checkpoint, but deal.II "
                             "did not detect its presence when you called `cmake'."));
#endif

    surface_pressure                = prm.get_double ("Surface pressure");
    adiabatic_surface_temperature   = prm.get_double ("Adiabatic surface temperature");
    pressure_normalization          = prm.get("Pressure normalization");

    use_operator_splitting          = prm.get_bool("Use operator splitting");

    prm.enter_subsection ("Mesh refinement");
    {
      initial_global_refinement    = prm.get_integer ("Initial global refinement");
      initial_adaptive_refinement  = prm.get_integer ("Initial adaptive refinement");

      adaptive_refinement_interval = prm.get_integer ("Time steps between mesh refinement");
      refinement_fraction          = prm.get_double ("Refinement fraction");
      coarsening_fraction          = prm.get_double ("Coarsening fraction");
      adapt_by_fraction_of_cells   = prm.get_bool ("Adapt by fraction of cells");
      min_grid_level               = prm.get_integer ("Minimum refinement level");

      AssertThrow(refinement_fraction >= 0 && coarsening_fraction >= 0,
                  ExcMessage("Refinement/coarsening fractions must be positive."));
      AssertThrow(refinement_fraction+coarsening_fraction <= 1,
                  ExcMessage("Refinement and coarsening fractions must be <= 1."));
      AssertThrow(min_grid_level <= initial_global_refinement,
                  ExcMessage("Minimum refinement level must not be larger than "
                             "Initial global refinement."));

      // extract the list of times at which additional refinement is requested
      // then sort it and convert it to seconds
      additional_refinement_times
        = Utilities::string_to_double
          (Utilities::split_string_list(prm.get ("Additional refinement times")));
      std::sort (additional_refinement_times.begin(),
                 additional_refinement_times.end());
      if (convert_to_years == true)
        for (unsigned int i=0; i<additional_refinement_times.size(); ++i)
          additional_refinement_times[i] *= year_in_seconds;

      skip_solvers_on_initial_refinement = prm.get_bool("Skip solvers on initial refinement");
      skip_setup_initial_conditions_on_initial_refinement = prm.get_bool("Skip setup initial conditions on initial refinement");

      if (skip_setup_initial_conditions_on_initial_refinement == true && skip_solvers_on_initial_refinement == false)
        AssertThrow(false, ExcMessage("Cannot execute solvers if no initial conditions are set up. "
                                      "You must set skip_solvers_on_initial_refinement to true."));

      run_postprocessors_on_initial_refinement = prm.get_bool("Run postprocessors on initial refinement");

      if (skip_setup_initial_conditions_on_initial_refinement == true && run_postprocessors_on_initial_refinement == true)
        AssertThrow(false, ExcMessage("Cannot run postprocessors if no initial conditions are set up. "
                                      "You must set run_postprocessors_on_initial_refinement to false."));
    }
    prm.leave_subsection ();

    prm.enter_subsection ("Postprocess");
    {
      run_postprocessors_on_nonlinear_iterations = prm.get_bool("Run postprocessors on nonlinear iterations");
    }
    prm.leave_subsection ();

    prm.enter_subsection ("Formulation");
    {
      // The following options each have a set of conditions to be met in order
      // for the formulation to be consistent, however, most of
      // the information is not available at this point. Therefore, the error checking is done
      // in Simulator<dim>::check_consistency_of_formulation() after the initialization of
      // material models, heating plugins, and adiabatic conditions.
      formulation = Formulation::parse(prm.get("Formulation"));
      if (formulation == Formulation::isentropic_compression)
        {
          formulation_mass_conservation = Formulation::MassConservation::isentropic_compression;
          formulation_temperature_equation = Formulation::TemperatureEquation::real_density;
        }
      else if (formulation == Formulation::boussinesq_approximation)
        {
          formulation_mass_conservation = Formulation::MassConservation::incompressible;
          formulation_temperature_equation = Formulation::TemperatureEquation::reference_density_profile;
        }
      else if (formulation == Formulation::anelastic_liquid_approximation)
        {
          // equally possible: implicit_reference_profile
          formulation_mass_conservation = Formulation::MassConservation::reference_density_profile;
          formulation_temperature_equation = Formulation::TemperatureEquation::reference_density_profile;
        }
      else if (formulation == Formulation::custom)
        {
          formulation_mass_conservation = Formulation::MassConservation::parse(prm.get("Mass conservation"));
          formulation_temperature_equation = Formulation::TemperatureEquation::parse(prm.get("Temperature equation"));
        }
      else AssertThrow(false, ExcNotImplemented());

      enable_additional_stokes_rhs = prm.get_bool ("Enable additional Stokes RHS");
      enable_elasticity = prm.get_bool("Enable elasticity");
      enable_prescribed_dilation = prm.get_bool("Enable prescribed dilation");
    }
    prm.leave_subsection ();


    prm.enter_subsection ("Melt settings");
    {
      include_melt_transport = prm.get_bool ("Include melt transport");
    }
    prm.leave_subsection();

    prm.enter_subsection ("Nullspace removal");
    {
      nullspace_removal = NullspaceRemoval::none;
      std::vector<std::string> nullspace_names =
        Utilities::split_string_list(prm.get("Remove nullspace"));
      AssertThrow(Utilities::has_unique_entries(nullspace_names),
                  ExcMessage("The list of strings for the parameter "
                             "'Nullspace removal/Remove nullspace' contains entries more than once. "
                             "This is not allowed. Please check your parameter file."));

      for (unsigned int i=0; i<nullspace_names.size(); ++i)
        {
          if (nullspace_names[i]=="net rotation")
            nullspace_removal = typename NullspaceRemoval::Kind(
                                  nullspace_removal | NullspaceRemoval::net_rotation);
          else if (nullspace_names[i]=="angular momentum")
            nullspace_removal = typename NullspaceRemoval::Kind(
                                  nullspace_removal | NullspaceRemoval::angular_momentum);
          else if (nullspace_names[i]=="net translation")
            nullspace_removal = typename NullspaceRemoval::Kind(
                                  nullspace_removal | NullspaceRemoval::net_translation_x |
                                  NullspaceRemoval::net_translation_y | ( dim == 3 ?
                                                                          NullspaceRemoval::net_translation_z : 0) );
          else if (nullspace_names[i]=="net x translation")
            nullspace_removal = typename NullspaceRemoval::Kind(
                                  nullspace_removal | NullspaceRemoval::net_translation_x);
          else if (nullspace_names[i]=="net y translation")
            nullspace_removal = typename NullspaceRemoval::Kind(
                                  nullspace_removal | NullspaceRemoval::net_translation_y);
          else if (nullspace_names[i]=="net z translation")
            nullspace_removal = typename NullspaceRemoval::Kind(
                                  nullspace_removal | NullspaceRemoval::net_translation_z);
          else if (nullspace_names[i]=="linear x momentum")
            nullspace_removal = typename NullspaceRemoval::Kind(
                                  nullspace_removal | NullspaceRemoval::linear_momentum_x);
          else if (nullspace_names[i]=="linear y momentum")
            nullspace_removal = typename NullspaceRemoval::Kind(
                                  nullspace_removal | NullspaceRemoval::linear_momentum_y);
          else if (nullspace_names[i]=="linear z momentum")
            nullspace_removal = typename NullspaceRemoval::Kind(
                                  nullspace_removal | NullspaceRemoval::linear_momentum_z);
          else if (nullspace_names[i]=="linear momentum")
            nullspace_removal = typename NullspaceRemoval::Kind(
                                  nullspace_removal | NullspaceRemoval::linear_momentum_x |
                                  NullspaceRemoval::linear_momentum_y | ( dim == 3 ?
                                                                          NullspaceRemoval::linear_momentum_z : 0) );
          else
            AssertThrow(false, ExcInternalError());
        }
    }
    prm.leave_subsection ();

    prm.enter_subsection ("Checkpointing");
    {
      checkpoint_time_secs = prm.get_integer ("Time between checkpoint");
      checkpoint_steps     = prm.get_integer ("Steps between checkpoint");

#ifndef DEAL_II_WITH_ZLIB
      AssertThrow ((checkpoint_time_secs == 0)
                   &&
                   (checkpoint_steps == 0),
                   ExcMessage ("You need to have deal.II configured with the `libz' "
                               "option if you want to generate checkpoints, but deal.II "
                               "did not detect its presence when you called `cmake'."));
#endif
    }
    prm.leave_subsection ();

    prm.enter_subsection ("Compositional fields");
    {
      n_compositional_fields = prm.get_integer ("Number of fields");
    }
    prm.leave_subsection();

    prm.enter_subsection ("Discretization");
    {
      stokes_velocity_degree = prm.get_integer ("Stokes velocity polynomial degree");
      temperature_degree     = prm.get_integer ("Temperature polynomial degree");
      composition_degree     = prm.get_integer ("Composition polynomial degree");
      use_locally_conservative_discretization
        = prm.get_bool ("Use locally conservative discretization");
      use_equal_order_interpolation_for_stokes
        = prm.get_bool ("Use equal order interpolation for Stokes");
      use_discontinuous_temperature_discretization
        = prm.get_bool("Use discontinuous temperature discretization");
      use_discontinuous_composition_discretization
        = prm.get_bool("Use discontinuous composition discretization");

      AssertThrow(use_discontinuous_composition_discretization == true || composition_degree > 0,
                  ExcMessage("Using a composition polynomial degree of 0 (cell-wise constant composition) "
                             "is only supported if a discontinuous composition discretization is selected."));

      prm.enter_subsection ("Stabilization parameters");
      {
        advection_stabilization_method = AdvectionStabilizationMethod::parse(prm.get("Stabilization method"));
        use_artificial_viscosity_smoothing  = prm.get_bool ("Use artificial viscosity smoothing");
        stabilization_alpha                 = prm.get_integer ("alpha");

        stabilization_c_R                   = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("cR"))),
                                                                                      n_compositional_fields+1,
                                                                                      "cR");
        stabilization_beta                  = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("beta"))),
                                                                                      n_compositional_fields+1,
                                                                                      "beta");

        stabilization_gamma                 = prm.get_double ("gamma");
        discontinuous_penalty               = prm.get_double ("Discontinuous penalty");
        use_limiter_for_discontinuous_temperature_solution
          = prm.get_bool("Use limiter for discontinuous temperature solution");
        use_limiter_for_discontinuous_composition_solution
          = prm.get_bool("Use limiter for discontinuous composition solution");
        global_temperature_max_preset       = prm.get_double ("Global temperature maximum");
        global_temperature_min_preset       = prm.get_double ("Global temperature minimum");
        global_composition_max_preset       = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double
                                                                                      (Utilities::split_string_list(prm.get ("Global composition maximum"))),
                                                                                      n_compositional_fields,
                                                                                      "Global composition maximum");
        global_composition_min_preset       = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double
                                                                                      (Utilities::split_string_list(prm.get ("Global composition minimum"))),
                                                                                      n_compositional_fields,
                                                                                      "Global composition minimum");
      }
      prm.leave_subsection ();

      AssertThrow ((use_locally_conservative_discretization ||
                    use_equal_order_interpolation_for_stokes)
                   ||
                   (stokes_velocity_degree > 1),
                   ExcMessage ("The polynomial degree for the velocity field "
                               "specified in the 'Stokes velocity polynomial degree' "
                               "parameter must be at least 2, unless you are using "
                               "the 'Use equal order interpolation for Stokes' parameter, "
                               "or if you are using "
                               "a locally conservative discretization as specified by the "
                               "'Use locally conservative discretization' parameter. "
                               "In both of these cases, the polynomial degree used for the "
                               "velocity may be equal to one."
                               "\n\n"
                               "The restriction exists because by default, the pressure element "
                               "is of one degree lower and continuous, and if you selected "
                               "a linear element for the velocity, you'd need a continuous "
                               "element of degree zero for the pressure, which does not exist. "
                               "On the other hand, if using equal-order interpolation, "
                               "choosing the polynomial degree as one yields a Q1-Q1 "
                               "element; using a locally conservative discretization "
                               "with polynomial degree of one yields a Q1-P0 "
                               "element."));

      AssertThrow (! (use_locally_conservative_discretization &&
                      use_equal_order_interpolation_for_stokes),
                   ExcMessage ("You have tried to use both the 'Use locally "
                               "conservative discretization' and 'Use equal order "
                               "interpolation for Stokes' parameters in the input "
                               "file. However, their use is incompatible: you "
                               "can only select one of the two."));
    }
    prm.leave_subsection ();

    prm.enter_subsection ("Temperature field");
    {
      std::string x_temperature_method
        = prm.get ("Temperature method");

      if (x_temperature_method == "field")
        temperature_method = AdvectionFieldMethod::fem_field;
      else if (x_temperature_method == "prescribed field")
        temperature_method = AdvectionFieldMethod::prescribed_field;
      else
        AssertThrow(false,ExcNotImplemented());
    }
    prm.leave_subsection();

    prm.enter_subsection ("Compositional fields");
    {
      if (include_melt_transport && (n_compositional_fields == 0))
        {
          AssertThrow (false,
                       ExcMessage ("If melt transport is included in the model, "
                                   "there has to be at least one compositional field."));
        }

      names_of_compositional_fields = Utilities::split_string_list (prm.get("Names of fields"));
      AssertThrow ((names_of_compositional_fields.size() == 0) ||
                   (names_of_compositional_fields.size() == n_compositional_fields),
                   ExcMessage ("The length of the list of names for the compositional "
                               "fields needs to either be empty or have length equal to "
                               "the number of compositional fields."));

      // check that the names use only allowed characters, are not empty strings and are unique
      for (unsigned int i=0; i<names_of_compositional_fields.size(); ++i)
        {
          AssertThrow (names_of_compositional_fields[i].find_first_not_of("abcdefghijklmnopqrstuvwxyz"
                                                                          "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                                                          "0123456789_") == std::string::npos,
                       ExcMessage("Invalid character in field " + names_of_compositional_fields[i] + ". "
                                  "Names of compositional fields should consist of a "
                                  "combination of letters, numbers and underscores."));
          AssertThrow (names_of_compositional_fields[i].size() > 0,
                       ExcMessage("Invalid name of field " + names_of_compositional_fields[i] + ". "
                                  "Names of compositional fields need to be non-empty."));
          for (unsigned int j=0; j<i; ++j)
            AssertThrow (names_of_compositional_fields[i] != names_of_compositional_fields[j],
                         ExcMessage("Names of compositional fields have to be unique! " + names_of_compositional_fields[i] +
                                    " is used more than once."));
        }

      // default names if list is empty
      if (names_of_compositional_fields.size() == 0)
        for (unsigned int i=0; i<n_compositional_fields; ++i)
          names_of_compositional_fields.push_back("C_" + Utilities::int_to_string(i+1));

      // if we want to solve the melt transport equations, check that one of the fields
      // has the name porosity
      if (include_melt_transport && std::find(names_of_compositional_fields.begin(),
                                              names_of_compositional_fields.end(), "porosity")
          == names_of_compositional_fields.end())
        {
          AssertThrow (false, ExcMessage ("If melt transport is included in the model, "
                                          "there has to be at least one compositional field "
                                          "with the name `porosity'."));
        }

      const std::vector<int> n_normalized_fields = Utilities::string_to_int
                                                   (Utilities::split_string_list(prm.get ("List of normalized fields")));
      normalized_fields = std::vector<unsigned int> (n_normalized_fields.begin(),
                                                     n_normalized_fields.end());

      AssertThrow (normalized_fields.size() <= n_compositional_fields,
                   ExcMessage("Invalid input parameter file: Too many entries in List of normalized fields"));

      std::vector<std::string> x_compositional_field_methods
        = Utilities::split_string_list
          (prm.get ("Compositional field methods"));

      AssertThrow ((x_compositional_field_methods.size() == 0) ||
                   (x_compositional_field_methods.size() == 1) ||
                   (x_compositional_field_methods.size() == n_compositional_fields),
                   ExcMessage ("The length of the list of names for the field method of compositional "
                               "fields needs to be empty, or have one entry, or have a length equal to "
                               "the number of compositional fields."));

      // If no method is specified set the default, which is solve every composition
      // by a continuous field method
      if (x_compositional_field_methods.size() == 0)
        x_compositional_field_methods = std::vector<std::string> (n_compositional_fields,"field");
      // If only one method is specified assume all fields are solved by this method
      else if (x_compositional_field_methods.size() == 1)
        x_compositional_field_methods = std::vector<std::string> (n_compositional_fields,x_compositional_field_methods[0]);


      // Parse all field methods and store them, the vector should be empty
      // since nobody should have written into it yet.
      Assert(compositional_field_methods.size() == 0,
             ExcInternalError());
      compositional_field_methods.resize(n_compositional_fields);
      for (unsigned int i = 0; i < n_compositional_fields; ++i)
        {
          if (x_compositional_field_methods[i] == "field")
            compositional_field_methods[i] = AdvectionFieldMethod::fem_field;
          else if (x_compositional_field_methods[i] == "particles")
            compositional_field_methods[i] = AdvectionFieldMethod::particles;
          else if (x_compositional_field_methods[i] == "volume of fluid")
            {
              AssertThrow (dim==2,
                           ExcMessage ("The 'volume of fluid' method is currently "
                                       "only implemented for two-dimensional "
                                       "computations."));
              compositional_field_methods[i] = AdvectionFieldMethod::volume_of_fluid;
            }
          else if (x_compositional_field_methods[i] == "static")
            compositional_field_methods[i] = AdvectionFieldMethod::static_field;
          else if (x_compositional_field_methods[i] == "melt field")
            compositional_field_methods[i] = AdvectionFieldMethod::fem_melt_field;
          else if (x_compositional_field_methods[i] == "prescribed field")
            compositional_field_methods[i] = AdvectionFieldMethod::prescribed_field;
          else if (x_compositional_field_methods[i] == "prescribed field with diffusion")
            compositional_field_methods[i] = AdvectionFieldMethod::prescribed_field_with_diffusion;
          else
            AssertThrow(false,ExcNotImplemented());
        }

      // Enable Volume of Fluid field tracking if any compositional_field_methods are volume_of_fluid
      volume_of_fluid_tracking_enabled =
        (std::count(compositional_field_methods.begin(),compositional_field_methods.end(),AdvectionFieldMethod::volume_of_fluid)
         > 0);

      if (std::find(compositional_field_methods.begin(), compositional_field_methods.end(), AdvectionFieldMethod::fem_melt_field)
          != compositional_field_methods.end())
        AssertThrow (this->include_melt_transport,
                     ExcMessage ("The advection method 'melt field' can only be selected if melt "
                                 "transport is used in the simulation."));

      const std::vector<std::string> x_mapped_particle_properties
        = Utilities::split_string_list
          (prm.get ("Mapped particle properties"));

      const unsigned int number_of_particle_fields =
        std::count(compositional_field_methods.begin(),compositional_field_methods.end(),AdvectionFieldMethod::particles);

      AssertThrow ((x_mapped_particle_properties.size() == number_of_particle_fields)
                   || (x_mapped_particle_properties.size() == 0),
                   ExcMessage ("The list of names for the mapped particle property fields needs to either be empty or have a length equal to "
                               "the number of compositional fields that are interpolated from particle properties."));

      for (const auto &p : x_mapped_particle_properties)
        {
          // each entry has the format (white space is optional):
          // <name> : <value (might have spaces)> [component]
          //
          // first tease apart the two halves
          const std::vector<std::string> split_parts = Utilities::split_string_list (p, ':');
          AssertThrow (split_parts.size() == 2,
                       ExcMessage ("The format for mapped particle properties  "
                                   "requires that each entry has the form `"
                                   "<name of field> : <particle property> [component]', "
                                   "but there does not appear to be a colon in the entry <"
                                   + p
                                   + ">."));

          // the easy part: get the name of the compositional field
          const std::string key = split_parts[0];

          // check that the names used are actually names of fields,
          // are solved by particles, and are unique in this list
          std::vector<std::string>::iterator field_name_iterator = std::find(names_of_compositional_fields.begin(),
                                                                             names_of_compositional_fields.end(), key);
          AssertThrow (field_name_iterator
                       != names_of_compositional_fields.end(),
                       ExcMessage ("Name of field <" + key +
                                   "> appears in the parameter "
                                   "<Compositional fields/Mapped particle properties>, but "
                                   "there is no field with this name."));

          const unsigned int compositional_field_index = std::distance(names_of_compositional_fields.begin(),
                                                                       field_name_iterator);

          AssertThrow (compositional_field_methods[compositional_field_index]
                       == AdvectionFieldMethod::particles,
                       ExcMessage ("The field <" + key +
                                   "> appears in the parameter <Compositional fields/Mapped particle properties>, but "
                                   "is not advected by a particle method."));

          AssertThrow (std::count(names_of_compositional_fields.begin(),
                                  names_of_compositional_fields.end(), key) == 1,
                       ExcMessage ("Name of field <" + key +
                                   "> appears more than once in the parameter "
                                   "<Compositional fields/Mapped particle properties>."));

          // now for the rest. since we don't know whether there is a
          // component selector, start reading at the end and subtract
          // a number that might be a component selector
          std::string particle_property = split_parts[1];
          std::string component;
          if ((particle_property.size()>3) &&
              (particle_property[particle_property.size()-1] == ']'))
            {
              particle_property.erase (--particle_property.end());

              // this handles the (rare) case of multi digit components
              while ((particle_property[particle_property.size()-1] >= '0') &&
                     (particle_property[particle_property.size()-1] <= '9'))
                {
                  component.insert(component.begin(),particle_property[particle_property.size()-1]);
                  particle_property.erase (--particle_property.end());
                }

              AssertThrow (particle_property[particle_property.size()-1] == '[',
                           ExcMessage("Problem in parsing a component selector from the string <"
                                      + split_parts[1] + ">. A component selector has to be of the "
                                      "form [x], where x must be an unsigned integer between 0 "
                                      "and the maximum number of components of this particle property."));

              particle_property.erase (--particle_property.end());
            }

          // we've stopped reading component selectors now.
          // eat spaces that may be at the end of particle_property to get key
          while ((particle_property.size()>0) && (particle_property[particle_property.size()-1] == ' '))
            particle_property.erase (--particle_property.end());

          // finally, put it into the list
          mapped_particle_properties.insert(std::make_pair(compositional_field_index,
                                                           std::make_pair(particle_property,atoi(component.c_str()))));
        }
    }
    prm.leave_subsection ();


    // now also get the parameter related to material model averaging
    prm.enter_subsection ("Material model");
    {
      material_averaging
        = MaterialModel::MaterialAveraging::parse_averaging_operation_name
          (prm.get ("Material averaging"));
    }
    prm.leave_subsection ();

    // then, finally, let user additions that do not go through the usual
    // plugin mechanism, declare their parameters if they have subscribed
    // to the relevant signals
    SimulatorSignals<dim>::parse_additional_parameters (*this, prm);

    AssertThrow((!use_direct_stokes_solver) || (nullspace_removal == NullspaceRemoval::none),
                ExcMessage("Because of the difference in system partitioning, nullspace removal is "
                           "currently not compatible with the direct solver. "
                           "Please turn off one or both of the options 'Nullspace removal/Remove nullspace', "
                           "or 'Use direct solver for Stokes system', or contribute code to enable "
                           "this feature combination."));
  }



  template <int dim>
  void
  Parameters<dim>::
  parse_geometry_dependent_parameters(ParameterHandler &prm,
                                      const GeometryModel::Interface<dim> &geometry_model)
  {
    prm.enter_subsection ("Mesh deformation");
    {
      // Test here for whether there are any boundary indicators active
      // for which mesh deformation objects are to be set.
      const std::vector<std::string> x_mesh_deformation_boundary_indicators
        = Utilities::split_string_list(prm.get("Mesh deformation boundary indicators"),";");
      mesh_deformation_enabled = !x_mesh_deformation_boundary_indicators.empty();
    }
    prm.leave_subsection();

    prm.enter_subsection ("Boundary traction model");
    {
      const std::vector<std::string> x_prescribed_traction_boundary_indicators
        = Utilities::split_string_list
          (prm.get ("Prescribed traction boundary indicators"));
      for (const auto &p : x_prescribed_traction_boundary_indicators)
        {
          // each entry has the format (white space is optional):
          // <id> [x][y][z] : <value (might have spaces)>
          //
          // first tease apart the two halves
          const std::vector<std::string> split_parts = Utilities::split_string_list (p, ':');
          AssertThrow (split_parts.size() == 2,
                       ExcMessage ("The format for prescribed traction boundary indicators "
                                   "requires that each entry has the form `"
                                   "<id> [x][y][z] : <value>', but there does not "
                                   "appear to be a colon in the entry <"
                                   + p
                                   + ">."));

          // the easy part: get the value
          const std::string value = split_parts[1];

          // now for the rest. since we don't know whether there is a
          // component selector, start reading at the end and subtracting
          // letters x, y and z
          std::string key_and_comp = split_parts[0];
          std::string comp;
          while ((key_and_comp.size()>0) &&
                 ((key_and_comp[key_and_comp.size()-1] == 'x')
                  ||
                  (key_and_comp[key_and_comp.size()-1] == 'y')
                  ||
                  ((key_and_comp[key_and_comp.size()-1] == 'z') && (dim==3))))
            {
              comp += key_and_comp[key_and_comp.size()-1];
              key_and_comp.erase (--key_and_comp.end());
            }

          // we've stopped reading component selectors now. there are three
          // possibilities:
          // - no characters are left. this means that key_and_comp only
          //   consisted of a single word that only consisted of 'x', 'y'
          //   and 'z's. then this would have been a mistake to classify
          //   as a component selector, and we better undo it
          // - the last character of key_and_comp is not a whitespace. this
          //   means that the last word in key_and_comp ended in an 'x', 'y'
          //   or 'z', but this was not meant to be a component selector.
          //   in that case, put these characters back.
          // - otherwise, we split successfully. eat spaces that may be at
          //   the end of key_and_comp to get key
          if (key_and_comp.size() == 0)
            key_and_comp.swap (comp);
          else if (key_and_comp[key_and_comp.size()-1] != ' ')
            {
              key_and_comp += comp;
              comp = "";
            }
          else
            {
              while ((key_and_comp.size()>0) && (key_and_comp[key_and_comp.size()-1] == ' '))
                key_and_comp.erase (--key_and_comp.end());
            }

          // finally, try to translate the key into a boundary_id. then
          // make sure we haven't seen it yet
          types::boundary_id boundary_id;
          try
            {
              boundary_id = geometry_model.translate_symbolic_boundary_name_to_id(key_and_comp);
            }
          catch (const std::string &error)
            {
              AssertThrow (false, ExcMessage ("While parsing the entry <Boundary traction model/Prescribed "
                                              "traction indicators>, there was an error. Specifically, "
                                              "the conversion function complained as follows: "
                                              + error));
            }

          AssertThrow (prescribed_traction_boundary_indicators.find(boundary_id)
                       == prescribed_traction_boundary_indicators.end(),
                       ExcMessage ("Boundary indicator <" + Utilities::int_to_string(boundary_id) +
                                   "> appears more than once in the list of indicators "
                                   "for nonzero traction boundaries."));

          // finally, put it into the list
          prescribed_traction_boundary_indicators[boundary_id] =
            std::pair<std::string,std::string>(comp,value);
        }
    }
    prm.leave_subsection ();

    prm.enter_subsection ("Boundary heat flux model");
    {
      try
        {
          const std::vector<types::boundary_id> x_fixed_heat_flux_boundary_indicators
            = geometry_model.translate_symbolic_boundary_names_to_ids(Utilities::split_string_list
                                                                      (prm.get ("Fixed heat flux boundary indicators")));
          fixed_heat_flux_boundary_indicators
            = std::set<types::boundary_id> (x_fixed_heat_flux_boundary_indicators.begin(),
                                            x_fixed_heat_flux_boundary_indicators.end());
        }
      catch (const std::string &error)
        {
          AssertThrow (false, ExcMessage ("While parsing the entry <Boundary heat flux model/Fixed heat flux "
                                          "boundary indicators>, there was an error. Specifically, "
                                          "the conversion function complained as follows: "
                                          + error));
        }
    }
    prm.leave_subsection ();
  }



  template <int dim>
  void Simulator<dim>::declare_parameters (ParameterHandler &prm)
  {
    Parameters<dim>::declare_parameters (prm);
    Melt::Parameters<dim>::declare_parameters (prm);
    Newton::Parameters::declare_parameters (prm);
    MeshDeformation::MeshDeformationHandler<dim>::declare_parameters (prm);
    Postprocess::Manager<dim>::declare_parameters (prm);
    MeshRefinement::Manager<dim>::declare_parameters (prm);
    TimeStepping::Manager<dim>::declare_parameters (prm);
    MaterialModel::declare_parameters<dim> (prm);
    HeatingModel::Manager<dim>::declare_parameters (prm);
    GeometryModel::declare_parameters <dim>(prm);
    InitialTopographyModel::declare_parameters <dim>(prm);
    GravityModel::declare_parameters<dim> (prm);
    InitialTemperature::Manager<dim>::declare_parameters (prm);
    InitialComposition::Manager<dim>::declare_parameters (prm);
    PrescribedStokesSolution::declare_parameters<dim> (prm);
    BoundaryTemperature::Manager<dim>::declare_parameters (prm);
    BoundaryComposition::Manager<dim>::declare_parameters (prm);
    AdiabaticConditions::declare_parameters<dim> (prm);
    BoundaryVelocity::Manager<dim>::declare_parameters (prm);
    BoundaryTraction::declare_parameters<dim> (prm);
    BoundaryHeatFlux::declare_parameters<dim> (prm);
  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template Parameters<dim>::Parameters (ParameterHandler &prm, \
                                        MPI_Comm mpi_communicator); \
  template void Parameters<dim>::declare_parameters (ParameterHandler &prm); \
  template void Parameters<dim>::parse_parameters(ParameterHandler &prm, \
                                                  const MPI_Comm mpi_communicator); \
  template void Parameters<dim>::parse_geometry_dependent_parameters(ParameterHandler &prm, \
                                                                     const GeometryModel::Interface<dim> &geometry_model); \
  template void Simulator<dim>::declare_parameters (ParameterHandler &prm);

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
