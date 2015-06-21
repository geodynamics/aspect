/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/simulator.h>
#include <aspect/global.h>

#include <deal.II/base/parameter_handler.h>

#include <dirent.h>
#include <stdlib.h>


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
                       Patterns::Integer (2,4),
                       "The number of space dimensions you want to run this program in. "
                       "ASPECT can run in 2 and 3 space dimensions.");
    prm.declare_entry ("Additional shared libraries", "",
                       Patterns::List (Patterns::FileName()),
                       "A list of names of additional shared libraries that should be loaded "
                       "upon starting up the program. The names of these files can contain absolute "
                       "or relative paths (relative to the directory in which you call ASPECT). "
                       "In fact, file names that are do not contain any directory "
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
                       Patterns::Bool (),
                       "A flag indicating whether the computation should be resumed from "
                       "a previously saved state (if true) or start from scratch (if false).");

#ifndef DEAL_II_WITH_ZLIB
    AssertThrow (resume_computation == false,
                 ExcMessage ("You need to have deal.II configured with the 'libz' "
                             "option if you want to resume a computation from a checkpoint, but deal.II "
                             "did not detect its presence when you called 'cmake'."));
#endif

    prm.declare_entry ("Max nonlinear iterations", "10",
                       Patterns::Integer (0),
                       "The maximal number of nonlinear iterations to be performed.");

    prm.declare_entry ("Start time", "0",
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
                       "If this flag is set to 'true' conversion to years happens; if "
                       "it is 'false', no such conversion happens. Note that when 'true', "
                       "some input such as prescribed velocities should also use years "
                       "instead of seconds.");

    prm.declare_entry ("CFL number", "1.0",
                       Patterns::Double (0),
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
                       Patterns::Double (0),
                       "Set a maximum time step size for the solver to use. Generally the time step "
                       "based on the CFL number should be sufficient, but for complicated models "
                       "or benchmarking it may be useful to limit the time step to some value. "
                       "The default value is a value so that when converted from years into seconds "
                       "it equals the largest number representable by a floating "
                       "point number, implying an unlimited time step."
                       "Units: Years or seconds, depending on the ``Use years "
                       "in output instead of seconds'' parameter.");

    prm.declare_entry ("Use conduction timestep", "false",
                       Patterns::Bool (),
                       "Mantle convection simulations are often focused on convection "
                       "dominated systems. However, these codes can also be used to "
                       "investigate systems where heat conduction plays a dominant role. "
                       "This parameter indicates whether the simulator should also use "
                       "heat conduction in determining the length of each time step.");

    prm.declare_entry ("Nonlinear solver scheme", "IMPES",
                       Patterns::Selection ("IMPES|iterated IMPES|iterated Stokes|Stokes only|Advection only"),
                       "The kind of scheme used to resolve the nonlinearity in the system. "
                       "'IMPES' is the classical IMplicit Pressure Explicit Saturation scheme "
                       "in which ones solves the temperatures and Stokes equations exactly "
                       "once per time step, one after the other. The 'iterated IMPES' scheme "
                       "iterates this decoupled approach by alternating the solution of the "
                       "temperature and Stokes systems. The 'iterated Stokes' scheme solves "
                       "the temperature equation once at the beginning of each time step "
                       "and then iterates out the solution of the Stokes equation. The 'Stokes only' "
                       "scheme only solves the Stokes system and ignores compositions and the "
                       "temperature equation (careful, the material model must not depend on "
                       "the temperature; mostly useful for Stokes benchmarks). The 'Advection only'"
                       "scheme only solves the temperature and other advection systems and instead "
                       "of solving for the Stokes system, a prescribed velocity and pressure is "
                       "used");

    prm.declare_entry ("Nonlinear solver tolerance", "1e-5",
                       Patterns::Double(0,1),
                       "A relative tolerance up to which the nonlinear solver "
                       "will iterate. This parameter is only relevant if "
                       "Nonlinear solver scheme is set to 'iterated Stokes' or "
                       "'iterated IMPES'.");

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
                       "of the domain is zero; the surface of the domain is determined by asking "
                       "the geometry model whether a particular face of the geometry has a zero "
                       "or small `depth'. If the value of this parameter is `volume' then the "
                       "pressure is normalized so that the domain average is zero. If `no' is "
                       "given, the no pressure normalization is performed.");

    prm.declare_entry ("Surface pressure", "0",
                       Patterns::Double(),
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

    prm.declare_entry ("Adiabatic surface temperature", "0",
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

    prm.declare_entry ("Use direct solver for Stokes system", "false",
                       Patterns::Bool(),
                       "If set to true the linear system for the Stokes equation will "
                       "be solved using Trilinos klu, otherwise an iterative Schur "
                       "complement solver is used. The direct solver is only efficient "
                       "for small problems.");

    prm.declare_entry ("Linear solver tolerance", "1e-7",
                       Patterns::Double(0,1),
                       "A relative tolerance up to which the linear Stokes systems in each "
                       "time or nonlinear step should be solved. The absolute tolerance will "
                       "then be $\\| M x_0 - F \\| \\cdot \\text{tol}$, where $x_0 = (0,p_0)$ "
                       "is the initial guess of the pressure, $M$ is the system matrix, "
                       "F is the right-hand side, and tol is the parameter specified here. "
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
    prm.declare_entry ("Number of cheap Stokes solver steps", "30",
                       Patterns::Integer(0),
                       "As explained in the ASPECT paper (Kronbichler, Heister, and Bangerth, "
                       "GJI 2012) we first try to solve the Stokes system in every time "
                       "step using a GMRES iteration with a poor but cheap "
                       "preconditioner. By default, we try whether we can converge the GMRES "
                       "solver in 30 such iterations before deciding that we need a better "
                       "preconditioner. This is sufficient for simple problems with constant "
                       "viscosity and we never need the second phase with the more expensive "
                       "preconditioner. On the other hand, for more complex problems, and in "
                       "particular for problems with strongly varying viscosity, the 30 "
                       "cheap iterations don't actually do very much good and one might skip "
                       "this part right away. In that case, this parameter can be set to "
                       "zero, i.e., we immediately start with the better but more expensive "
                       "preconditioner.");

    prm.declare_entry ("Temperature solver tolerance", "1e-12",
                       Patterns::Double(0,1),
                       "The relative tolerance up to which the linear system for "
                       "the temperature system gets solved. See 'linear solver "
                       "tolerance' for more details.");

    prm.declare_entry ("Composition solver tolerance", "1e-12",
                       Patterns::Double(0,1),
                       "The relative tolerance up to which the linear system for "
                       "the composition system gets solved. See 'linear solver "
                       "tolerance' for more details.");

    // next declare parameters that pertain to the equations to be
    // solved, along with boundary conditions etc. note that at this
    // point we do not know yet which geometry model we will use, so
    // we do not know which symbolic names will be valid to address individual
    // parts of the boundary. we can only work around this by allowing any string
    // to indicate a boundary
    prm.enter_subsection ("Model settings");
    {
      prm.declare_entry ("Fixed temperature boundary indicators", "",
                         Patterns::List (Patterns::Anything()),
                         "A comma separated list of names denoting those boundaries "
                         "on which the temperature is fixed and described by the "
                         "boundary temperature object selected in its own section "
                         "of this input file. All boundary indicators used by the geometry "
                         "but not explicitly listed here will end up with no-flux "
                         "(insulating) boundary conditions."
                         "\n\n"
                         "The names of the boundaries listed here can either by "
                         "numeric numbers (in which case they correspond to the numerical "
                         "boundary indicators assigned by the geometry object), or they "
                         "can correspond to any of the symbolic names the geometry object "
                         "may have provided for each part of the boundary. You may want "
                         "to compare this with the documentation of the geometry model you "
                         "use in your model."
                         "\n\n"
                         "This parameter only describes which boundaries have a fixed "
                         "temperature, but not what temperature should hold on these "
                         "boundaries. The latter piece of information needs to be "
                         "implemented in a plugin in the BoundaryTemperature "
                         "group, unless an existing implementation in this group "
                         "already provides what you want.");
      prm.declare_entry ("Fixed composition boundary indicators", "",
                         Patterns::List (Patterns::Anything()),
                         "A comma separated list of names denoting those boundaries "
                         "on which the composition is fixed and described by the "
                         "boundary composition object selected in its own section "
                         "of this input file. All boundary indicators used by the geometry "
                         "but not explicitly listed here will end up with no-flux "
                         "(insulating) boundary conditions."
                         "\n\n"
                         "The names of the boundaries listed here can either by "
                         "numeric numbers (in which case they correspond to the numerical "
                         "boundary indicators assigned by the geometry object), or they "
                         "can correspond to any of the symbolic names the geometry object "
                         "may have provided for each part of the boundary. You may want "
                         "to compare this with the documentation of the geometry model you "
                         "use in your model."
                         "\n\n"
                         "This parameter only describes which boundaries have a fixed "
                         "composition, but not what composition should hold on these "
                         "boundaries. The latter piece of information needs to be "
                         "implemented in a plugin in the BoundaryComposition "
                         "group, unless an existing implementation in this group "
                         "already provides what you want.");
      prm.declare_entry ("Zero velocity boundary indicators", "",
                         Patterns::List (Patterns::Anything()),
                         "A comma separated list of names denoting those boundaries "
                         "on which the velocity is zero."
                         "\n\n"
                         "The names of the boundaries listed here can either by "
                         "numeric numbers (in which case they correspond to the numerical "
                         "boundary indicators assigned by the geometry object), or they "
                         "can correspond to any of the symbolic names the geometry object "
                         "may have provided for each part of the boundary. You may want "
                         "to compare this with the documentation of the geometry model you "
                         "use in your model.");
      prm.declare_entry ("Tangential velocity boundary indicators", "",
                         Patterns::List (Patterns::Anything()),
                         "A comma separated list of names denoting those boundaries "
                         "on which the velocity is tangential and unrestrained, i.e., free-slip where "
                         "no external forces act to prescribe a particular tangential "
                         "velocity (although there is a force that requires the flow to "
                         "be tangential)."
                         "\n\n"
                         "The names of the boundaries listed here can either by "
                         "numeric numbers (in which case they correspond to the numerical "
                         "boundary indicators assigned by the geometry object), or they "
                         "can correspond to any of the symbolic names the geometry object "
                         "may have provided for each part of the boundary. You may want "
                         "to compare this with the documentation of the geometry model you "
                         "use in your model.");
      prm.declare_entry ("Free surface boundary indicators", "",
                         Patterns::List (Patterns::Anything()),
                         "A comma separated list of names denoting those boundaries "
                         "where there is a free surface. Set to nothing to disable all "
                         "free surface computations."
                         "\n\n"
                         "The names of the boundaries listed here can either by "
                         "numeric numbers (in which case they correspond to the numerical "
                         "boundary indicators assigned by the geometry object), or they "
                         "can correspond to any of the symbolic names the geometry object "
                         "may have provided for each part of the boundary. You may want "
                         "to compare this with the documentation of the geometry model you "
                         "use in your model.");
      prm.declare_entry ("Prescribed velocity boundary indicators", "",
                         Patterns::Map (Patterns::Anything(),
                                        Patterns::Selection(VelocityBoundaryConditions::get_names<dim>())),
                         "A comma separated list denoting those boundaries "
                         "on which the velocity is prescribed, i.e., where unknown "
                         "external forces act to prescribe a particular velocity. This is "
                         "often used to prescribe a velocity that equals that of "
                         "overlying plates."
                         "\n\n"
                         "The format of valid entries for this parameter is that of a map "
                         "given as ``key1 [selector]: value1, key2 [selector]: value2, key3: value3, ...'' where "
                         "each key must be a valid boundary indicator (which is either an "
                         "integer or the symbolic name the geometry model in use may have "
                         "provided for this part of the boundary) "
                         "and each value must be one of the currently implemented boundary "
                         "velocity models. ``selector'' is an optional string given as a subset "
                         "of the letters 'xyz' that allows you to apply the boundary conditions "
                         "only to the components listed. As an example, '1 y: function' applies "
                         "the type 'function' to the y component on boundary 1. Without a selector "
                         "it will affect all components of the velocity."
                         "\n\n"
                         "Note that the no-slip boundary condition is "
                         "a special case of the current one where the prescribed velocity "
                         "happens to be zero. It can thus be implemented by indicating that "
                         "a particular boundary is part of the ones selected "
                         "using the current parameter and using ``zero velocity'' as "
                         "the boundary values. Alternatively, you can simply list the "
                         "part of the boundary on which the velocity is to be zero with "
                         "the parameter ``Zero velocity boundary indicator'' in the "
                         "current parameter section."
                         "\n\n"
                         "Note that when ``Use years in output instead of seconds'' is set "
                         "to true, velocity should be given in m/yr. ");
      prm.declare_entry ("Prescribed traction boundary indicators", "",
                         Patterns::Map (Patterns::Anything(),
                                        Patterns::Selection(TractionBoundaryConditions::get_names<dim>())),
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
                         "of the letters 'xyz' that allows you to apply the boundary conditions "
                         "only to the components listed. As an example, '1 y: function' applies "
                         "the type 'function' to the y component on boundary 1. Without a selector "
                         "it will affect all components of the traction.");
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
    prm.leave_subsection ();

    prm.enter_subsection ("Mesh refinement");
    {
      prm.declare_entry ("Initial global refinement", "2",
                         Patterns::Integer (0),
                         "The number of global refinement steps performed on "
                         "the initial coarse mesh, before the problem is first "
                         "solved there.");
      prm.declare_entry ("Initial adaptive refinement", "2",
                         Patterns::Integer (0),
                         "The number of adaptive refinement steps performed after "
                         "initial global refinement but while still within the first "
                         "time step.");
      prm.declare_entry ("Time steps between mesh refinement", "10",
                         Patterns::Integer (0),
                         "The number of time steps after which the mesh is to be "
                         "adapted again based on computed error indicators. If 0 "
                         "then the mesh will never be changed.");
      prm.declare_entry ("Refinement fraction", "0.3",
                         Patterns::Double(0,1),
                         "The fraction of cells with the largest error that "
                         "should be flagged for refinement.");
      prm.declare_entry ("Coarsening fraction", "0.05",
                         Patterns::Double(0,1),
                         "The fraction of cells with the smallest error that "
                         "should be flagged for coarsening.");
      prm.declare_entry ("Minimum refinement level", "0",
                         Patterns::Integer (0),
                         "The minimum refinement level each cell should have, "
                         "and that can not be exceeded by coarsening. "
                         "Should not be higher than the 'Initial global refinement' "
                         "parameter.");
      prm.declare_entry ("Additional refinement times", "",
                         Patterns::List (Patterns::Double(0)),
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
                         "Whether or not the postproccessors should be run at the end "
                         "of each of ths initial adaptive refinement cycles at the "
                         "of the simulation start.");
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
                         "pair conform with the usual LBB (Babuska-Brezzi) condition. In "
                         "other words, we are using a Taylor-Hood element for the Stoeks "
                         "equations and this parameter indicates the polynomial degree of it. "
                         "Units: None.");
      prm.declare_entry ("Temperature polynomial degree", "2",
                         Patterns::Integer (1),
                         "The polynomial degree to use for the temperature variable. "
                         "Units: None.");
      prm.declare_entry ("Composition polynomial degree", "2",
                         Patterns::Integer (1),
                         "The polynomial degree to use for the composition variable(s). "
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
                         "polynomial space $P_ {-q}$ of polynomials of degree $q$ in "
                         "each variable separately. Here, $q$ is one less than the value "
                         "given in the parameter ``Stokes velocity polynomial degree''. "
                         "As a consequence of choosing this "
                         "element, it can be shown if the medium is considered incompressible "
                         "that the computed discrete velocity "
                         "field $\\mathbf u_h$ satisfies the property $\\int_ {\\partial K} \\mathbf u_h "
                         "\\cdot \\mathbf n = 0$ for every cell $K$, i.e., for each cell inflow and "
                         "outflow exactly balance each other as one would expect for an "
                         "incompressible medium. In other words, the velocity field is locally "
                         "conservative."
                         "\n\n"
                         "On the other hand, if this parameter is "
                         "set to ``false'', then the finite element space is chosen as $Q_q$. "
                         "This choice does not yield the local conservation property but "
                         "has the advantage of requiring fewer degrees of freedom. Furthermore, "
                         "the error is generally smaller with this choice."
                         "\n\n"
                         "For an in-depth discussion of these issues and a quantitative evaluation "
                         "of the different choices, see \\cite {KHB12} .");

      prm.enter_subsection ("Stabilization parameters");
      {
        prm.declare_entry ("alpha", "2",
                           Patterns::Integer (1, 2),
                           "The exponent $\\alpha$ in the entropy viscosity stabilization. Valid "
                           "options are 1 or 2. The recommended setting is 2. (This parameter does "
                           "not correspond to any variable in the 2012 GJI paper by Kronbichler, "
                           "Heister and Bangerth that describes ASPECT. Rather, the paper always uses "
                           "2 as the exponent in the definition of the entropy, following eq. (15).)."
                           "Units: None.");
        prm.declare_entry ("cR", "0.33",
                           Patterns::Double (0),
                           "The $c_R$ factor in the entropy viscosity "
                           "stabilization. (For historical reasons, the name used here is different "
                           "from the one used in the 2012 GJI paper by Kronbichler, "
                           "Heister and Bangerth that describes ASPECT. This parameter corresponds "
                           "to the factor $\\alpha_E$ in the formulas following equation (15) of "
                           "the paper. After further experiments, we have also chosen to use a "
                           "different value than described there.) Units: None.");
        prm.declare_entry ("beta", "0.078",
                           Patterns::Double (0),
                           "The $\\beta$ factor in the artificial viscosity "
                           "stabilization. An appropriate value for 2d is 0.078 "
                           "and 0.117 for 3d. (For historical reasons, the name used here is different "
                           "from the one used in the 2012 GJI paper by Kronbichler, "
                           "Heister and Bangerth that describes ASPECT. This parameter corresponds "
                           "to the factor $\\alpha_\\text {max}$ in the formulas following equation (15) of "
                           "the paper. After further experiments, we have also chosen to use a "
                           "different value than described there: It can be chosen as stated there for "
                           "uniformly refined meshes, but it needs to be chosen larger if the mesh has "
                           "cells that are not squares or cubes.) Units: None.");
      }
      prm.leave_subsection ();
    }
    prm.leave_subsection ();

    prm.enter_subsection ("Compositional fields");
    {
      prm.declare_entry ("Number of fields", "0",
                         Patterns::Integer (0),
                         "The number of fields that will be advected along with the flow field, excluding "
                         "velocity, pressure and temperature.");
      prm.declare_entry ("Names of fields", "",
                         Patterns::List(Patterns::Anything()),
                         "A user-defined name for each of the compositional fields requested.");
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
                         "The process of averaging, and where it may be used, is "
                         "discussed in more detail in "
                         "Section~\\ref{sec:sinker-with-averaging}."
                         "\n\n"
                         "More averaging schemes are available in the averaging material "
                         "model. This material model is a ``compositing material model'' "
                         "which can be used in combination with other material models.");
    }
    prm.leave_subsection ();

    // also declare the parameters that the FreeSurfaceHandler needs
    Simulator<dim>::FreeSurfaceHandler::declare_parameters (prm);

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

    resume_computation      = prm.get_bool ("Resume computation");
    CFL_number              = prm.get_double ("CFL number");
    use_conduction_timestep = prm.get_bool ("Use conduction timestep");
    convert_to_years        = prm.get_bool ("Use years in output instead of seconds");
    timing_output_frequency = prm.get_integer ("Timing output frequency");

    maximum_time_step       = prm.get_double("Maximum time step");
    if (convert_to_years == true)
      maximum_time_step *= year_in_seconds;


    if (prm.get ("Nonlinear solver scheme") == "IMPES")
      nonlinear_solver = NonlinearSolver::IMPES;
    else if (prm.get ("Nonlinear solver scheme") == "iterated IMPES")
      nonlinear_solver = NonlinearSolver::iterated_IMPES;
    else if (prm.get ("Nonlinear solver scheme") == "iterated Stokes")
      nonlinear_solver = NonlinearSolver::iterated_Stokes;
    else if (prm.get ("Nonlinear solver scheme") == "Stokes only")
      nonlinear_solver = NonlinearSolver::Stokes_only;
    else if (prm.get ("Nonlinear solver scheme") == "Advection only")
      nonlinear_solver = NonlinearSolver::Advection_only;
    else
      AssertThrow (false, ExcNotImplemented());

    nonlinear_tolerance = prm.get_double("Nonlinear solver tolerance");

    max_nonlinear_iterations = prm.get_integer ("Max nonlinear iterations");
    start_time              = prm.get_double ("Start time");
    if (convert_to_years == true)
      start_time *= year_in_seconds;

    output_directory        = prm.get ("Output directory");
    if (output_directory.size() == 0)
      output_directory = "./";
    else if (output_directory[output_directory.size()-1] != '/')
      output_directory += "/";

    // verify that the output directory actually exists. if it doesn't, create
    // it on processor zero
    if ((Utilities::MPI::this_mpi_process(mpi_communicator) == 0) &&
        (opendir(output_directory.c_str()) == NULL))
      {
        std::cout << "\n"
                  << "-----------------------------------------------------------------------------\n"
                  << "The output directory <" << output_directory
                  << "> provided in the input file appears not to exist.\n"
                  << "ASPECT will create it for you.\n"
                  << "-----------------------------------------------------------------------------\n\n"
                  << std::endl;

        // create the directory. we could call the 'mkdir()' function directly, but
        // this can only create a single level of directories. if someone has specified
        // a nested subdirectory as output directory, and if multiple parts of the path
        // do not exist, this would fail. working around this is easiest by just calling
        // 'mkdir -p' from the command line
        const int error = system ((std::string("mkdir -p '") + output_directory + "'").c_str());

        AssertThrow (error==0,
                     ExcMessage (std::string("Can't create the output directory at <") + output_directory + ">"));
      }

    surface_pressure              = prm.get_double ("Surface pressure");
    adiabatic_surface_temperature = prm.get_double ("Adiabatic surface temperature");
    pressure_normalization        = prm.get("Pressure normalization");

    use_direct_stokes_solver      = prm.get_bool("Use direct solver for Stokes system");
    linear_stokes_solver_tolerance= prm.get_double ("Linear solver tolerance");
    n_cheap_stokes_solver_steps   = prm.get_integer ("Number of cheap Stokes solver steps");
    temperature_solver_tolerance  = prm.get_double ("Temperature solver tolerance");
    composition_solver_tolerance  = prm.get_double ("Composition solver tolerance");

    prm.enter_subsection ("Mesh refinement");
    {
      initial_global_refinement   = prm.get_integer ("Initial global refinement");
      initial_adaptive_refinement = prm.get_integer ("Initial adaptive refinement");

      adaptive_refinement_interval= prm.get_integer ("Time steps between mesh refinement");
      refinement_fraction         = prm.get_double ("Refinement fraction");
      coarsening_fraction         = prm.get_double ("Coarsening fraction");
      min_grid_level              = prm.get_integer ("Minimum refinement level");

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

      run_postprocessors_on_initial_refinement = prm.get_bool("Run postprocessors on initial refinement");
    }
    prm.leave_subsection ();

    prm.enter_subsection ("Model settings");
    {
      {
        nullspace_removal = NullspaceRemoval::none;
        std::vector<std::string> nullspace_names =
          Utilities::split_string_list(prm.get("Remove nullspace"));
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
              nullspace_removal = typename       NullspaceRemoval::Kind(
                                    nullspace_removal | NullspaceRemoval::linear_momentum_x);
            else if (nullspace_names[i]=="linear y momentum")
              nullspace_removal = typename       NullspaceRemoval::Kind(
                                    nullspace_removal | NullspaceRemoval::linear_momentum_y);
            else if (nullspace_names[i]=="linear z momentum")
              nullspace_removal = typename       NullspaceRemoval::Kind(
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
                   ExcMessage ("You need to have deal.II configured with the 'libz' "
                               "option if you want to generate checkpoints, but deal.II "
                               "did not detect its presence when you called 'cmake'."));
#endif
    }
    prm.leave_subsection ();

    prm.enter_subsection ("Discretization");
    {
      stokes_velocity_degree = prm.get_integer ("Stokes velocity polynomial degree");
      temperature_degree     = prm.get_integer ("Temperature polynomial degree");
      composition_degree     = prm.get_integer ("Composition polynomial degree");
      use_locally_conservative_discretization
        = prm.get_bool ("Use locally conservative discretization");

      prm.enter_subsection ("Stabilization parameters");
      {
        stabilization_alpha = prm.get_integer ("alpha");
        stabilization_c_R   = prm.get_double ("cR");
        stabilization_beta  = prm.get_double ("beta");
      }
      prm.leave_subsection ();

      AssertThrow (use_locally_conservative_discretization ||
                   (stokes_velocity_degree > 1),
                   ExcMessage ("The polynomial degree for the velocity field "
                               "specified in the 'Stokes velocity polynomial degree' "
                               "parameter must be at least 2, unless you are using "
                               "a locally conservative discretization as specified by the "
                               "'Use locally conservative discretization' parameter. "
                               "This is because in the former case, the pressure element "
                               "is of one degree lower and continuous, and if you selected "
                               "a linear element for the velocity, you'd need a continuous "
                               "element of degree zero for the pressure, which does not exist."))
    }
    prm.leave_subsection ();

    prm.enter_subsection ("Compositional fields");
    {
      n_compositional_fields = prm.get_integer ("Number of fields");

      names_of_compositional_fields = Utilities::split_string_list (prm.get("Names of fields"));
      AssertThrow ((names_of_compositional_fields.size() == 0) ||
                   (names_of_compositional_fields.size() == n_compositional_fields),
                   ExcMessage ("The length of the list of names for the compositional "
                               "fields needs to either be empty or have length equal to "
                               "the number of compositional fields."));

      // check that the names use only allowed characters, are not empty strings and are unique
      for (unsigned int i=0; i<names_of_compositional_fields.size(); ++i)
        {
          Assert (names_of_compositional_fields[i].find_first_not_of("abcdefghijklmnopqrstuvwxyz"
                                                                     "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                                                     "0123456789_") == std::string::npos,
                  ExcMessage("Invalid character in field " + names_of_compositional_fields[i] + ". "
                             "Names of compositional fields should consist of a "
                             "combination of letters, numbers and underscores."));
          Assert (names_of_compositional_fields[i].size() > 0,
                  ExcMessage("Invalid name of field " + names_of_compositional_fields[i] + ". "
                             "Names of compositional fields need to be non-empty."));
          for (unsigned int j=0; j<i; ++j)
            Assert (names_of_compositional_fields[i] != names_of_compositional_fields[j],
                    ExcMessage("Names of compositional fields have to be unique! " + names_of_compositional_fields[i] +
                               " is used more than once."));
        }

      // default names if list is empty
      if (names_of_compositional_fields.size() == 0)
        for (unsigned int i=0; i<n_compositional_fields; ++i)
          names_of_compositional_fields.push_back("C_" + Utilities::int_to_string(i+1));

      const std::vector<int> n_normalized_fields = Utilities::string_to_int
                                                   (Utilities::split_string_list(prm.get ("List of normalized fields")));
      normalized_fields = std::vector<unsigned int> (n_normalized_fields.begin(),
                                                     n_normalized_fields.end());

      AssertThrow (normalized_fields.size() <= n_compositional_fields,
                   ExcMessage("Invalid input parameter file: Too many entries in List of normalized fields"));
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
  }



  template <int dim>
  void
  Parameters<dim>::
  parse_geometry_dependent_parameters(ParameterHandler &prm,
                                      const GeometryModel::Interface<dim> &geometry_model)
  {
    prm.enter_subsection ("Model settings");
    {
      try
        {
          const std::vector<types::boundary_id> x_fixed_temperature_boundary_indicators
            = geometry_model.translate_symbolic_boundary_names_to_ids(Utilities::split_string_list
                                                                      (prm.get ("Fixed temperature boundary indicators")));
          fixed_temperature_boundary_indicators
            = std::set<types::boundary_id> (x_fixed_temperature_boundary_indicators.begin(),
                                            x_fixed_temperature_boundary_indicators.end());
        }
      catch (const std::string &error)
        {
          AssertThrow (false, ExcMessage ("While parsing the entry <Model settings/Fixed temperature "
                                          "boundary indicators>, there was an error. Specifically, "
                                          "the conversion function complained as follows: "
                                          + error));
        }

      try
        {
          const std::vector<types::boundary_id> x_fixed_composition_boundary_indicators
            = geometry_model.translate_symbolic_boundary_names_to_ids (Utilities::split_string_list
                                                                       (prm.get ("Fixed composition boundary indicators")));
          fixed_composition_boundary_indicators
            = std::set<types::boundary_id> (x_fixed_composition_boundary_indicators.begin(),
                                            x_fixed_composition_boundary_indicators.end());
        }
      catch (const std::string &error)
        {
          AssertThrow (false, ExcMessage ("While parsing the entry <Model settings/Fixed composition "
                                          "boundary indicators>, there was an error. Specifically, "
                                          "the conversion function complained as follows: "
                                          + error));
        }

      try
        {
          const std::vector<types::boundary_id> x_zero_velocity_boundary_indicators
            = geometry_model.translate_symbolic_boundary_names_to_ids(Utilities::split_string_list
                                                                      (prm.get ("Zero velocity boundary indicators")));
          zero_velocity_boundary_indicators
            = std::set<types::boundary_id> (x_zero_velocity_boundary_indicators.begin(),
                                            x_zero_velocity_boundary_indicators.end());
        }
      catch (const std::string &error)
        {
          AssertThrow (false, ExcMessage ("While parsing the entry <Model settings/Zero velocity "
                                          "boundary indicators>, there was an error. Specifically, "
                                          "the conversion function complained as follows: "
                                          + error));
        }

      try
        {
          const std::vector<types::boundary_id> x_tangential_velocity_boundary_indicators
            = geometry_model.translate_symbolic_boundary_names_to_ids(Utilities::split_string_list
                                                                      (prm.get ("Tangential velocity boundary indicators")));
          tangential_velocity_boundary_indicators
            = std::set<types::boundary_id> (x_tangential_velocity_boundary_indicators.begin(),
                                            x_tangential_velocity_boundary_indicators.end());
        }
      catch (const std::string &error)
        {
          AssertThrow (false, ExcMessage ("While parsing the entry <Model settings/Tangential velocity "
                                          "boundary indicators>, there was an error. Specifically, "
                                          "the conversion function complained as follows: "
                                          + error));
        }

      try
        {
          const std::vector<types::boundary_id> x_free_surface_boundary_indicators
            = geometry_model.translate_symbolic_boundary_names_to_ids(Utilities::split_string_list
                                                                      (prm.get ("Free surface boundary indicators")));
          free_surface_boundary_indicators
            = std::set<types::boundary_id> (x_free_surface_boundary_indicators.begin(),
                                            x_free_surface_boundary_indicators.end());

          free_surface_enabled = !free_surface_boundary_indicators.empty();
        }
      catch (const std::string &error)
        {
          AssertThrow (false, ExcMessage ("While parsing the entry <Model settings/Free surface "
                                          "boundary indicators>, there was an error. Specifically, "
                                          "the conversion function complained as follows: "
                                          + error));
        }

      const std::vector<std::string> x_prescribed_velocity_boundary_indicators
        = Utilities::split_string_list
          (prm.get ("Prescribed velocity boundary indicators"));
      for (std::vector<std::string>::const_iterator p = x_prescribed_velocity_boundary_indicators.begin();
           p != x_prescribed_velocity_boundary_indicators.end(); ++p)
        {
          // each entry has the format (white space is optional):
          // <id> [x][y][z] : <value (might have spaces)>
          //
          // first tease apart the two halves
          const std::vector<std::string> split_parts = Utilities::split_string_list (*p, ':');
          AssertThrow (split_parts.size() == 2,
                       ExcMessage ("The format for prescribed velocity boundary indicators "
                                   "requires that each entry has the form `"
                                   "<id> [x][y][z] : <value>', but there does not "
                                   "appear to be a colon in the entry <"
                                   + *p
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
              AssertThrow (false, ExcMessage ("While parsing the entry <Model settings/Prescribed "
                                              "velocity indicators>, there was an error. Specifically, "
                                              "the conversion function complained as follows: "
                                              + error));
            }

          AssertThrow (prescribed_velocity_boundary_indicators.find(boundary_id)
                       == prescribed_velocity_boundary_indicators.end(),
                       ExcMessage ("Boundary indicator <" + Utilities::int_to_string(boundary_id) +
                                   "> appears more than once in the list of indicators "
                                   "for nonzero velocity boundaries."));

          // finally, put it into the list
          prescribed_velocity_boundary_indicators[boundary_id] =
            std::pair<std::string,std::string>(comp,value);
        }

      const std::vector<std::string> x_prescribed_traction_boundary_indicators
        = Utilities::split_string_list
          (prm.get ("Prescribed traction boundary indicators"));
      for (std::vector<std::string>::const_iterator p = x_prescribed_traction_boundary_indicators.begin();
           p != x_prescribed_traction_boundary_indicators.end(); ++p)
        {
          // each entry has the format (white space is optional):
          // <id> [x][y][z] : <value (might have spaces)>
          //
          // first tease apart the two halves
          const std::vector<std::string> split_parts = Utilities::split_string_list (*p, ':');
          AssertThrow (split_parts.size() == 2,
                       ExcMessage ("The format for prescribed traction boundary indicators "
                                   "requires that each entry has the form `"
                                   "<id> [x][y][z] : <value>', but there does not "
                                   "appear to be a colon in the entry <"
                                   + *p
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
              AssertThrow (false, ExcMessage ("While parsing the entry <Model settings/Prescribed "
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
  }



  template <int dim>
  void Simulator<dim>::declare_parameters (ParameterHandler &prm)
  {
    Parameters<dim>::declare_parameters (prm);
    Postprocess::Manager<dim>::declare_parameters (prm);
    MeshRefinement::Manager<dim>::declare_parameters (prm);
    TerminationCriteria::Manager<dim>::declare_parameters (prm);
    MaterialModel::declare_parameters<dim> (prm);
    HeatingModel::Manager<dim>::declare_parameters (prm);
    GeometryModel::declare_parameters <dim>(prm);
    GravityModel::declare_parameters<dim> (prm);
    InitialConditions::declare_parameters<dim> (prm);
    CompositionalInitialConditions::declare_parameters<dim> (prm);
    PrescribedStokesSolution::declare_parameters<dim> (prm);
    BoundaryTemperature::declare_parameters<dim> (prm);
    BoundaryComposition::declare_parameters<dim> (prm);
    AdiabaticConditions::declare_parameters<dim> (prm);
    VelocityBoundaryConditions::declare_parameters<dim> (prm);
    TractionBoundaryConditions::declare_parameters<dim> (prm);
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
}
