/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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


#ifndef _aspect_simulator_h
#define _aspect_simulator_h

#include <deal.II/base/timer.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/symmetric_tensor.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/base/tensor_function.h>

DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <aspect/global.h>
#include <aspect/simulator_access.h>
#include <aspect/lateral_averaging.h>
#include <aspect/simulator_signals.h>
#include <aspect/material_model/interface.h>
#include <aspect/heating_model/interface.h>
#include <aspect/geometry_model/initial_topography_model/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/boundary_temperature/interface.h>
#include <aspect/boundary_heat_flux/interface.h>
#include <aspect/boundary_convective_heating/interface.h>
#include <aspect/boundary_composition/interface.h>
#include <aspect/initial_temperature/interface.h>
#include <aspect/initial_composition/interface.h>
#include <aspect/prescribed_stokes_solution/interface.h>
#include <aspect/boundary_velocity/interface.h>
#include <aspect/boundary_fluid_pressure/interface.h>
#include <aspect/boundary_traction/interface.h>
#include <aspect/mesh_refinement/interface.h>
#include <aspect/time_stepping/interface.h>
#include <aspect/postprocess/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/particle/manager.h>
#include <aspect/advection_field.h>

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>

#include <memory>
#include <thread>

namespace WorldBuilder
{
  class World;
}


namespace aspect
{
  template <int dim>
  class MeltHandler;

  template <int dim>
  class NewtonHandler;

  template <int dim>
  class StokesMatrixFreeHandler;

  namespace StokesSolver
  {
    template <int dim>
    class Direct;
  }

  template <int dim, int velocity_degree>
  class StokesMatrixFreeHandlerLocalSmoothingImplementation;

  namespace MeshDeformation
  {
    template <int dim>
    class MeshDeformationHandler;
  }

  template <int dim>
  class VolumeOfFluidHandler;

  namespace internal
  {
    namespace Assembly
    {
      namespace Scratch
      {
        template <int dim>      struct StokesPreconditioner;
        template <int dim>      struct StokesSystem;
        template <int dim>      struct AdvectionSystem;
      }

      namespace CopyData
      {
        template <int dim>      struct StokesPreconditioner;
        template <int dim>      struct StokesSystem;
        template <int dim>      struct AdvectionSystem;
      }
    }
  }

  namespace Assemblers
  {
    template <int dim>      class Interface;
    template <int dim>      class Manager;
  }

  struct DefectCorrectionResiduals
  {
    double initial_residual;
    double velocity_residual;
    double pressure_residual;
    double residual;
    double residual_old;
    double switch_initial_residual;
    double newton_residual_for_derivative_scaling_factor;
    std::pair<double,double> stokes_residuals;
  };

  /**
   * A data structure with all properties relevant to compute angular momentum and rotation.
   */
  template <int dim>
  struct RotationProperties
  {
    RotationProperties()
      :
      scalar_moment_of_inertia(numbers::signaling_nan<double>()),
      scalar_angular_momentum(numbers::signaling_nan<double>()),
      scalar_rotation(numbers::signaling_nan<double>()),
      tensor_moment_of_inertia(numbers::signaling_nan<SymmetricTensor<2,dim>>()),
      tensor_angular_momentum(numbers::signaling_nan<Tensor<1,dim>>()),
      tensor_rotation(numbers::signaling_nan<Tensor<1,dim>>())
    {};

    /**
     * Scalar properties for the two-dimensional case
     * with a fixed rotation axis (z).
     */
    double scalar_moment_of_inertia;
    double scalar_angular_momentum;
    double scalar_rotation;

    /**
     * Tensor properties for the three-dimensional case.
     */
    SymmetricTensor<2,dim> tensor_moment_of_inertia;
    Tensor<1,dim> tensor_angular_momentum;
    Tensor<1,dim> tensor_rotation;
  };

  /**
   * Exception to be thrown when the nonlinear solver needs too many iterations to converge.
   */
  DeclExceptionMsg(ExcNonlinearSolverNoConvergence,
                   "Nonlinear solver failed to converge in the prescribed number of steps. "
                   "Consider changing `Max nonlinear iterations` or `Nonlinear solver failure "
                   "strategy`.");

  /**
   * This is the main class of ASPECT. It implements the overall simulation
   * algorithm using the numerical methods discussed in the papers and manuals
   * that accompany ASPECT.
   *
   * @ingroup Simulator
   */
  template <int dim>
  class Simulator
  {
    public:
      /**
       * Constructor.
       *
       * @param mpi_communicator The MPI communicator on which this class is
       * to work. The class creates a clone of the actual communicator to make
       * its communications private from the rest of the world.
       *
       * @param prm The run-time parameter object from which this class
       * obtains its settings.
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      Simulator (const MPI_Comm mpi_communicator,
                 ParameterHandler &prm);

      /**
       * Destructor. Destroy what needs to be destroyed after waiting for all
       * threads that may still be doing something in the background.
       */
      ~Simulator ();

      /**
       * Declare the run-time parameters this class takes, and call the
       * respective <code>declare_parameters</code> functions of the
       * namespaces that describe geometries, material models, etc.
       *
       * @param prm The object in which the run-time parameters are to be
       * declared.
       *
       * This function is implemented in
       * <code>source/simulator/parameters.cc</code>.
       */
      static
      void declare_parameters (ParameterHandler &prm);

      /**
       * The function that runs the overall algorithm. It contains the loop
       * over all time steps as well as the logic of what to do when before
       * the loop starts and within the time loop.
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      void run ();

      /**
       * Write a connection graph of all of the plugins we know about, in the
       * format that the programs dot and neato understand. This allows for a
       * visualization of how all of the plugins that ASPECT knows about are
       * interconnected, and connect to other parts of the ASPECT code.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       *
       * @param output_stream The stream to write the output to.
       */
      void
      write_plugin_graph (std::ostream &output_stream) const;

      /**
       * Import Nonlinear Solver type.
       */
      using NonlinearSolver = typename Parameters<dim>::NonlinearSolver;

      /**
       * Import nullspace removal type.
       */
      using NullspaceRemoval = typename Parameters<dim>::NullspaceRemoval;

      using AdvectionField = aspect::AdvectionField;

    private:

      /**
       * A member variable that tracks whether we are completely done
       * with initialization and have started the time loop. This
       * variable needs to be the first one that is initialized in
       * the constructor, so that its value correctly tracks the
       * status of the overall object. As a consequence, it has to
       * be the first one declared in this class.
       *
       * The variable is set to @p true just before we start the time
       * stepping loop, but may be temporarily reset to @p false
       * if, for example, we are during the initial mesh refinement
       * steps where we start the time loop, but then go back to
       * initialization steps (mesh refinement, interpolation of initial
       * conditions, etc.) before re-starting the time loop.
       *
       * This variable is queried by
       * SimulatorAccess::simulator_is_past_initialization().
       */
      bool simulator_is_past_initialization;

      /**
       * A class that is empty but that can be used as a member variable and
       * whose constructor will be run in the order in which the member
       * variables are initialized. Because this class has a constructor that
       * takes a function object that it will execute whenever the member
       * variable is initialized, this allows running arbitrary actions in
       * between member variable initializers, for example if some member
       * variable is partially initialized at point A within the member
       * variable initializer list, its initialization can only be finalized
       * after point B (because it depends on what another member variable
       * decides to do), but needs to be finished by point C within the member
       * initialization. In such a case, one may have a member variable of the
       * current time placed in the list of member variables such that it is
       * initialized at point B, and then initialize it using a function
       * object that performs the finalization of initialization.
       */
      struct IntermediaryConstructorAction
      {
        IntermediaryConstructorAction (const std::function<void ()> &action);
      };

      /**
       * @name Top-level functions in the overall flow of the numerical
       * algorithm
       * @{
       */

      /**
       * The function that sets up the DoFHandler objects, It also sets up the
       * various partitioners and computes those constraints on the Stokes
       * variable and temperature that are the same between all time steps.
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      void setup_dofs ();

      /**
       * This function initializes the variables of the introspection object.
       * It is called by setup_dofs() right after distributing degrees of
       * freedom since this is when all of the information is available.
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      void setup_introspection ();

      /**
       * A function that is responsible for initializing the
       * temperature/compositional field before the first time step. The
       * temperature field then serves as the temperature from which the
       * velocity is computed during the first time step, and is subsequently
       * overwritten by the temperature field one gets by advancing by one
       * time step.
       *
       * This function is implemented in
       * <code>source/simulator/initial_conditions.cc</code>.
       */
      void set_initial_temperature_and_compositional_fields ();

      /**
       * A function that initializes the pressure variable before the first
       * time step. It does so by either interpolating (for continuous
       * pressure finite elements) or projecting (for discontinuous elements)
       * the adiabatic pressure computed from the material model.
       *
       * Note that the pressure so set is overwritten by the pressure in fact
       * computed during the first time step. We need this function, however,
       * so that the evaluation of pressure-dependent coefficients (e.g.
       * pressure dependent densities or thermal coefficients) during the
       * first time step has some useful pressure to start with.
       *
       * This function is implemented in
       * <code>source/simulator/initial_conditions.cc</code>.
       */
      void compute_initial_pressure_field ();

      /**
       * Fill the given @p constraints with constraints coming from the velocity boundary
       * conditions that do not change over time. This function is used by
       * setup_dofs();
       */
      void compute_initial_velocity_boundary_constraints (AffineConstraints<double> &constraints);

      /**
       * Fill the given @p constraints with constraints coming from the velocity boundary
       * conditions that do can change over time. This function is used by
       * compute_current_constraints().
       */
      void compute_current_velocity_boundary_constraints (AffineConstraints<double> &constraints);

      /**
       * Given the 'constraints' member that contains all constraints that are
       * independent of the time (e.g., hanging node constraints, tangential
       * flow constraints, etc), copy it over to 'current_constraints' and add
       * to the latter all constraints that do depend on time such as
       * temperature or velocity Dirichlet boundary conditions. This function
       * is therefore called at the beginning of every time step in
       * start_timestep(), but also when setting up the initial values.
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      void compute_current_constraints ();

      /**
       * Compute the factor by which we scale the second of
       * the Stokes equations (the "pressure scaling factor").
       * We compute the factor by taking the logarithmic
       * average of the viscosities we find on the cells in
       * this domain and dividing this average by a reference
       * length scale provided by the used geometry model.
       *
       * This function returns the pressure scaling variable.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      double compute_pressure_scaling_factor () const;

      /**
       * Do some housekeeping at the beginning of each time step. This
       * includes generating some screen output, adding some information to
       * the statistics file, and interpolating time-dependent boundary
       * conditions specific to this particular time step (the time
       * independent boundary conditions, for example for hanging nodes or for
       * tangential flow, are computed only once per mesh in setup_dofs()).
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      void start_timestep ();

      /**
       * Do the various steps necessary to assemble and solve the things
       * necessary in each time step.
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      void solve_timestep ();

      /**
       * This function implements one scheme for the various
       * steps necessary to assemble and solve the nonlinear problem.
       *
       * The `no Advection, no Stokes' scheme skips solving the temperature,
       * composition and Stokes equations, which permits to go directly to
       * postprocessing after setting up the initial condition.
       *
       * This function is implemented in
       * <code>source/simulator/solver_schemes.cc</code>.
       */
      void solve_no_advection_no_stokes ();

      /**
       * This function implements one scheme for the various
       * steps necessary to assemble and solve the nonlinear problem.
       *
       * The `no Advection, single Stokes' scheme only solves the Stokes system and
       * ignores compositions and the temperature equation.
       *
       * This function is implemented in
       * <code>source/simulator/solver_schemes.cc</code>.
       */
      void solve_no_advection_single_stokes ();

      /**
       * This function implements one scheme for the various
       * steps necessary to assemble and solve the nonlinear problem.
       *
       * The `first timestep only, single Stokes' scheme only solves the Stokes system,
       * for the initial timestep. This results in a `steady state' velocity field for
       * particle calculations.
       *
       * This function is implemented in
       * <code>source/simulator/solver_schemes.cc</code>.
       */
      void solve_no_advection_single_stokes_first_timestep_only ();

      /**
       * This function implements one scheme for the various
       * steps necessary to assemble and solve the nonlinear problem.
       *
       * The `no Advection, iterated Stokes' scheme only solves the Stokes system and
       * ignores compositions and the temperature equation (careful, the material
       * model must not depend on the temperature; mostly useful for
       * Stokes benchmarks).
       *
       * This function is implemented in
       * <code>source/simulator/solver_schemes.cc</code>.
       */
      void solve_no_advection_iterated_stokes ();

      /**
       * This function implements one scheme for the various
       * steps necessary to assemble and solve the nonlinear problem.
       *
       * The `no Advection, iterated defect correction Stokes' scheme
       * does not solve the temperature and composition equations
       * but only iterates out the solution of the Stokes
       * equation using Defect Correction (DC) Picard iterations.
       *
       * This function is implemented in
       * <code>source/simulator/solver_schemes.cc</code>.
       */
      void solve_no_advection_iterated_defect_correction_stokes ();

      /**
       * This function implements one scheme for the various
       * steps necessary to assemble and solve the nonlinear problem.
       *
       * The `single Advection, no Stokes' scheme only solves the temperature and other
       * advection systems and instead of solving for the Stokes system,
       * a prescribed velocity and pressure is used."
       *
       * This function is implemented in
       * <code>source/simulator/solver_schemes.cc</code>.
       */
      void solve_single_advection_no_stokes ();

      /**
       * This function implements one scheme for the various
       * steps necessary to assemble and solve the nonlinear problem.
       *
       * If `single Advection, single Stokes' is selected as the nonlinear solver scheme,
       * no nonlinear iterations are done, and the temperature, compositional fields and
       * Stokes equations are solved exactly once per time step, one after the other.
       *
       * This function is implemented in
       * <code>source/simulator/solver_schemes.cc</code>.
       */
      void solve_single_advection_single_stokes ();

      /**
       * This function implements one scheme for the various
       * steps necessary to assemble and solve the nonlinear problem.
       *
       * The `single Advection, iterated Stokes' scheme solves the temperature and
       * composition equations once at the beginning of each time step
       * and then iterates out the solution of the Stokes equation using
       * Picard iterations.
       *
       * This function is implemented in
       * <code>source/simulator/solver_schemes.cc</code>.
       */
      void solve_single_advection_iterated_stokes ();

      /**
       * This function implements one scheme for the various
       * steps necessary to assemble and solve the nonlinear problem.
       *
       * The `single Advection, iterated defect correction Stokes' scheme
       * solves the temperature and composition equations once at the beginning
       * of each time step and then iterates out the solution of the Stokes
       * equation using Defect Correction (DC) Picard iterations.
       *
       * This function is implemented in
       * <code>source/simulator/solver_schemes.cc</code>.
       */
      void solve_single_advection_iterated_defect_correction_stokes ();

      /**
       * This function implements one scheme for the various
       * steps necessary to assemble and solve the nonlinear problem.
       *
       * The `single Advection, iterated Newton Stokes' scheme solves the temperature and
       * composition equations once at the beginning of each time step
       * and then iterates out the solution of the Stokes equation using Newton iterations.
       * For the Stokes system it is able to switch from a defect correction form of
       * Picard iterations to Newton iterations after a certain tolerance or
       * number of iterations is reached. This can greatly improve the
       * convergence rate for particularly nonlinear viscosities.
       *
       * @param use_newton_iterations Sets whether this function should only use defect
       * correction iterations (use_newton_iterations = false) or also use Newton iterations
       * (use_newton_iterations = true).
       *
       * This function is implemented in
       * <code>source/simulator/solver_schemes.cc</code>.
       */
      void solve_single_advection_iterated_newton_stokes (bool use_newton_iterations);

      /**
       * This function implements one scheme for the various
       * steps necessary to assemble and solve the nonlinear problem.
       *
       * The `iterated Advection and Stokes' scheme iterates
       * by alternating the solution of the temperature, composition and Stokes systems.
       * This is essentially a type of Picard iterations for the whole
       * system of equations.
       *
       * This function is implemented in
       * <code>source/simulator/solver_schemes.cc</code>.
       */
      void solve_iterated_advection_and_stokes ();

      /**
       * This function implements one scheme for the various
       * steps necessary to assemble and solve the nonlinear problem.
       *
       * The `iterated Advection and defect correction Stokes' scheme
       * iterates over both the temperature and composition equations
       * and the Stokes equations at the same time. The Stokes equation
       * is iterated out using the Defect Correction (DC) form of the
       * Picard iterations.
       *
       * This function is implemented in
       * <code>source/simulator/solver_schemes.cc</code>.
       */
      void solve_iterated_advection_and_defect_correction_stokes ();

      /**
       * This function implements one scheme for the various
       * steps necessary to assemble and solve the nonlinear problem.
       *
       * The `iterated Advection and Newton Stokes' scheme iterates over solving the temperature,
       * composition, and Stokes equations just like `iterated Advection and Stokes', but
       * for the Stokes system it is able to switch from a defect correction form of
       * Picard iterations to Newton iterations after a certain tolerance or
       * number of iterations is reached. This can greatly improve the
       * convergence rate for particularly nonlinear viscosities.
       *
       * @param use_newton_iterations Sets whether this function should only use defect
       * correction iterations (use_newton_iterations = false) or also use Newton iterations
       * (use_newton_iterations = true).
       *
       * This function is implemented in
       * <code>source/simulator/solver_schemes.cc</code>.
       */
      void solve_iterated_advection_and_newton_stokes (bool use_newton_iterations);

      /**
       * Initiate the assembly of the Stokes preconditioner matrix via
       * assemble_stokes_preconditioner(), then set up the data structures to
       * actually build a preconditioner from this matrix.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void build_stokes_preconditioner ();

      /**
       * Initialize the preconditioner for the advection equation of field
       * index.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void build_advection_preconditioner (const AdvectionField &advection_field,
                                           aspect::LinearAlgebra::PreconditionILU &preconditioner,
                                           const double diagonal_strengthening);

      /**
       * Initiate the assembly of the Stokes matrix and right hand side.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void assemble_stokes_system ();

      /**
       * Assemble and solve the temperature equation.
       * This function returns the residual after solving.
       *
       * If the `residual` argument is not a `nullptr`, the function computes
       * the residual and puts it into this variable. The function returns
       * the current residual divided by the initial residual given as the
       * first argument. The two arguments may point to the same variable,
       * in which case the function first computes the residual and at
       * the end scales that residual by itself, thus returning 1.0.
       *
       * This function is implemented in
       * <code>source/simulator/solver_schemes.cc</code>.
       */
      double assemble_and_solve_temperature (const double &initial_residual = 0,
                                             double *residual = nullptr);

      /**
       * Solve the composition equations with whatever method is selected
       * (fields or particles). This function returns the residuals for
       * all fields after solving.
       *
       * If the `residual` argument is not a `nullptr`, the function computes
       * the residual and puts it into this variable. The function returns
       * the current residual divided by the initial residual given as the
       * first argument. The two arguments may point to the same variable,
       * in which case the function first computes the residual and at
       * the end scales that residual by itself, thus returning 1.0.
       *
       * This function is implemented in
       * <code>source/simulator/solver_schemes.cc</code>.
       */
      std::vector<double> assemble_and_solve_composition (const std::vector<double> &initial_residual = {},
                                                          std::vector<double> *residual = nullptr);

      /**
       * Assemble and solve the Stokes equation.
       * This function returns the nonlinear residual after solving.
       *
       * If the `residual` argument is not a `nullptr`, the function computes
       * the residual and puts it into this variable. The function returns
       * the current residual divided by the initial residual given as the
       * first argument. The two arguments may point to the same variable,
       * in which case the function first computes the residual and at
       * the end scales that residual by itself, thus returning 1.0.
       *
       * This function is implemented in
       * <code>source/simulator/solver_schemes.cc</code>.
       */
      double assemble_and_solve_stokes (const double &initial_nonlinear_residual = 0,
                                        double *nonlinear_residual = nullptr);

      /**
       * Do one step of the defect correction form of the Stokes equation;
       * i.e., assemble and solve the defect correction equations.
       * This function takes a structure of DefectCorrectionResiduals which
       * contains information about different residuals. The information in
       * this structure is updated by this function. The parameter use_picard
       * forces the use of the defect correction Picard iteration, if set to 'true' no
       * Newton derivatives are added to the matrix.
       *
       * This function is implemented in
       * <code>source/simulator/solver_schemes.cc</code>.
       */
      void do_one_defect_correction_Stokes_step(DefectCorrectionResiduals &dcr,
                                                const bool use_picard);

      /**
       * Initiate the assembly of one advection matrix and right hand side and
       * build a preconditioner for the matrix.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void assemble_advection_system (const AdvectionField &advection_field);

      /**
       * Solve one block of the temperature/composition linear system.
       * Return the initial nonlinear residual, i.e., if the linear system to
       * be solved is $Ax=b$, then return $\|Ax_0-b\|$ where $x_0$ is the
       * initial guess for the solution variable and is taken from the
       * current_linearization_point member variable.
       *
       * This function is implemented in
       * <code>source/simulator/solver.cc</code>.
       */
      double solve_advection (const AdvectionField &advection_field);

      /**
       * Interpolate the corresponding particle properties into the given
       * @p advection_fields solution fields.
       */
      void interpolate_particle_properties (const std::vector<AdvectionField> &advection_fields);

      /**
       * Solve the Stokes linear system.
       *
       * The function returns two pieces of information as a pair of doubles:
       * - The initial nonlinear residual, i.e., if the linear system to be
       *   solved is $Ax_{k+1}=b$, then we use $\|Ax_k-b\|$ where $x_k$ is the
       *   initial guess for the solution variable and is taken from
       *   the @p current_linearization_point member variable. For the
       *   purpose of this function, this residual is computed
       *   only from the velocity and pressure equations (i.e., for the 2x2 block
       *   system involving the velocity and pressure variables). A rationale
       *   for why this number is computed is given below.
       * - The final linear residual, i.e., if the linear system to be
       *   solved is $Ax_{k+1}=b$, then we use $\|Ax_{k+1}-b\|$ where $x_{k+1}$
       *   is the solution just computed. If we use a direct solver to compute
       *   the solution of the linear system, then this linear residual is of
       *   course zero (or at least quite close to it) and the function just
       *   returns a zero value without even attempting to compute the actual
       *   value. On the other hand, if the function uses an iterative solver,
       *   then the value of the final linear residual is related to the
       *   tolerance with which we solve the linear system and generally
       *   indicates how accurately or inaccurately the linear system has been
       *   solved.
       *
       * The two values are used in nonlinear solver schemes to assess how
       * accurate the solution was before the current solve (for the first
       * element of the returned pair) and how accurately the next iteration
       * will have to be solved (for the second element of the pair) when using
       * the Eisenstat-Walker method.
       *
       * @note If this function is called from a nonlinear solver -- e.g., the
       * `single Advection, iterated Stokes', or the
       * `iterated Advection and Stokes' solvers schemes --, then the
       * @p current_linearization_point is the solution of the previous
       * iteration (or the solution extrapolated from the previous time
       * steps, if this is the first nonlinear iteration). Let us call
       * this solution $x_k$ as above, where $x$ is a two-component
       * block vector that consists of velocity and pressure. This function
       * then assumes that we have already built the system matrix $A_k=A(x_k)$
       * and $F_k=F(x_k)$, both linearized around the previous solution. The
       * function solves the linear system $A_k x_{k+1} = F_k$ for the
       * solution $x_{k+1}$. If the linear system were solved exactly, then
       * that would imply that the <i>linear residual</i>
       * $\|A_k x_{k+1} - F_k\|$ were zero, or at least small. In other words,
       * its size does not tell us anything about how accurately we have
       * solved the nonlinear system. On the other hand, the <i>nonlinear
       * residual</i> $\|A_k x_k - F_k\|$ tells us something about how
       * accurately the previous guess $x_k$ already solved the nonlinear
       * system. Consequently, this is what this function returns. (In some
       * sense, this is not really what we are interested in: it tells us
       * how accurate the solution <i>already</i> was, and if it was already
       * pretty accurate, then we may not want to actually solve for
       * $x_{k+1}$. But, this would require that this function receives a
       * tolerance so that it can bail out early without actually solving
       * for $x_{k+1}$ if the tolerance is already reached. This function does
       * not actually do that -- in some sense, one may argue that if we have
       * already built the matrix and right hand side, we may as well solve
       * with them, whether or not the solution was already good. If it
       * happens to have been good already, then it will be even better after
       * the solve. If it was not good enough yet, then we have to solve
       * anyway.) In contrast to all of this, if we are using a Newton
       * solver, then $x_{k+1}$ is actually the Newton <i>update</i>
       * vector, for which we have no initial guess other than the zero
       * vector. In this case, the function simply returns $\|F_k\|$ as the
       * first element of the pair, where $F_k=F(x_k)$ is the residual
       * vector for the previous solution $x_k$.
       *
       * @param solution_vector The solution vector that is computed by this
       * function. This vector is a block vector that has the same block
       * structure as the full solution vector and its pressure and velocity
       * blocks will be overwritten by the solution of the Stokes system.
       *
       *
       * This function is implemented in
       * <code>source/simulator/solver.cc</code>.
       */
      std::pair<double,double>
      solve_stokes (LinearAlgebra::BlockVector &solution_vector);

      /**
       * This function is called at the end of every time step. It runs all
       * the postprocessors that have been listed in the input parameter file
       * (see the manual) in turn. In particular, this usually includes
       * generating graphical output every few time steps.
       *
       * The function also updates the statistics output file at the end of
       * each time step.
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      void postprocess ();

      /**
       * Refine the mesh according to error indicators calculated by
       * compute_refinement_criterion(), set up all necessary data structures
       * on this new mesh, and interpolate the old solutions onto the new
       * mesh.
       *
       * @param[in] max_grid_level The maximum refinement level of the mesh.
       * This is the sum of the initial global refinement and the initial
       * adaptive refinement (as provided by the user in the input file) and
       * in addition it gets increased by one at each additional refinement
       * time.
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      void refine_mesh (const unsigned int max_grid_level);

      /**
       * @}
       */

      /**
       * @name Functions used in saving the state of the program and
       * restarting from a saved state
       * @{
       */

      /**
       * Determine the id of the last good snapshot that was written by reading
       * the last_good_checkpoint.txt file from the output/checkpoint/ folder.
       * It will return numbers::invalid_unsigned_int if no snapshot exists.
       */
      unsigned int determine_last_good_snapshot() const;

      /**
       * Save the state of this program to a set of files in the output
       * directory. In reality, however, only some variables are stored (in
       * particular the mesh, the solution vectors, etc) whereas others can
       * either be re-generated (matrices, DoFHandler objects, etc) or are
       * read from the input parameter file. See the manual for more
       * information.
       *
       * This function is implemented in
       * <code>source/simulator/checkpoint_restart.cc</code>.
       */
      void create_snapshot();

      /**
       * Restore the state of this program from a set of files in the output
       * directory. In reality, however, only some variables are stored (in
       * particular the mesh, the solution vectors, etc) whereas others can
       * either be re-generated (matrices, DoFHandler objects, etc) or are
       * read from the input parameter file. See the manual for more
       * information. This function only restores those variables that can
       * neither be re-generated from other information nor are read from the
       * input parameter file.
       *
       * This function is implemented in
       * <code>source/simulator/checkpoint_restart.cc</code>.
       */
      void resume_from_snapshot();

      /**
       * Save a number of variables using BOOST serialization mechanism.
       *
       * This function is implemented in
       * <code>source/simulator/checkpoint_restart.cc</code>.
       */
      template <class Archive>
      void serialize (Archive &ar, const unsigned int version);
      /**
       * @}
       */

      /**
       * @name Functions used in setting up linear systems
       * @{
       */
      /**
       * Determine which of the components of our finite-element
       * system couple to each other. Depending on which equations
       * are solved and which solver is used this varies widely.
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      Table<2,DoFTools::Coupling>
      setup_system_matrix_coupling () const;

      /**
       * Set up the size and structure of the matrix used to store the
       * elements of the linear system.
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      void setup_system_matrix (const std::vector<IndexSet> &system_partitioning);

      /**
       * Set up the size and structure of the matrix used to store the
       * elements of the matrix that is used to build the
       * preconditioner for the system. This matrix is only used for
       * the Stokes system, so while it has the size of the whole
       * system, it only has entries in the velocity and pressure
       * blocks.
       *
       * This function is implemented in
       * <code>source/simulator/core.cc</code>.
       */
      void setup_system_preconditioner (const std::vector<IndexSet> &system_partitioning);

      /**
       * @}
       */

      /**
       * @name Functions, classes, and variables used in the assembly of linear systems
       * @{
       */

      /**
       * A member variable that stores, for the current simulation, what
       * functions need to be called in order to assemble linear systems,
       * matrices, and right hand side vectors.
       *
       * One would probably want this variable to just be a member of type
       * Assemblers::Manager<dim>, but this requires that
       * this type is declared in the current scope, and that would require
       * including <simulator/assemblers/interface.h> which we don't want because it's big.
       * Consequently, we just store a pointer to such an object, and create
       * the object pointed to at the top of set_assemblers().
       */
      std::unique_ptr<Assemblers::Manager<dim>> assemblers;

      /**
       * Determine, based on the run-time parameters of the current simulation,
       * which functions need to be called in order to assemble linear systems,
       * matrices, and right hand side vectors.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void set_assemblers ();

      /**
       * Determine, based on the run-time parameters of the current simulation,
       * which functions need to be called in order to assemble linear systems,
       * matrices, and right hand side vectors for the advection. This function
       * is used by both the default full Stokes solver and the Newton solvers,
       * but not by the two-phase flow solver.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void set_advection_assemblers ();

      /**
       * Determine, based on the run-time parameters of the current simulation,
       * which functions need to be called in order to assemble linear systems,
       * matrices, and right hand side vectors for the default full Stokes solver
       * i.e. without considering two-phase flow, or Newton solvers.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void set_stokes_assemblers ();

      /**
       * Initiate the assembly of the preconditioner for the Stokes system.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void assemble_stokes_preconditioner ();

      /**
       * Compute the integrals for the preconditioner for the Stokes system on
       * a single cell.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      local_assemble_stokes_preconditioner (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                            internal::Assembly::Scratch::StokesPreconditioner<dim> &scratch,
                                            internal::Assembly::CopyData::StokesPreconditioner<dim> &data);

      /**
       * Copy the contribution to the preconditioner for the Stokes system
       * from a single cell into the global matrix that stores these elements.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      copy_local_to_global_stokes_preconditioner (const internal::Assembly::CopyData::StokesPreconditioner<dim> &data);

      /**
       * Compute the integrals for the Stokes matrix and right hand side on a
       * single cell.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      local_assemble_stokes_system (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                    internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                    internal::Assembly::CopyData::StokesSystem<dim> &data);

      /**
       * Copy the contribution to the Stokes system from a single cell into
       * the global matrix that stores these elements.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      copy_local_to_global_stokes_system (const internal::Assembly::CopyData::StokesSystem<dim> &data);

      /**
       * Compute the integrals for one advection matrix and right hand side on
       * the faces of a single cell.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      local_assemble_advection_face_terms(const AdvectionField &advection_field,
                                          const typename DoFHandler<dim>::active_cell_iterator &cell,
                                          internal::Assembly::Scratch::AdvectionSystem<dim> &scratch,
                                          internal::Assembly::CopyData::AdvectionSystem<dim> &data);
      /**
       * Compute the integrals for one advection matrix and right hand side on
       * a single cell.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      local_assemble_advection_system (const AdvectionField &advection_field,
                                       const Vector<double>           &viscosity_per_cell,
                                       const typename DoFHandler<dim>::active_cell_iterator &cell,
                                       internal::Assembly::Scratch::AdvectionSystem<dim>  &scratch,
                                       internal::Assembly::CopyData::AdvectionSystem<dim> &data);

      /**
       * Copy the contribution to the advection system from a single cell into
       * the global matrix that stores these elements.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      copy_local_to_global_advection_system (const AdvectionField &advection_field,
                                             const internal::Assembly::CopyData::AdvectionSystem<dim> &data);

      /**
       * @}
       */

      /**
       * @name Helper functions
       * @{
       */
      /**
       * This routine adjusts the second block of the right hand side of a
       * Stokes system  (containing the term that comes from compressibility,
       * so that the system becomes compatible: $0=\int div u = \int g$. The
       * vector to adjust is given as the argument of this function. This
       * function makes use of the helper vector
       * pressure_shape_function_integrals that contains $h_i=(q_i,1)$ with
       * the pressure functions $q_i$ and we adjust the right hand side $g$ by
       * $h_i \int g / |\Omega|$.
       *
       * The purpose of this function is described in the second paper on the
       * numerical methods in ASPECT.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void make_pressure_rhs_compatible(LinearAlgebra::BlockVector &vector);

      /**
       * Fills a vector with the artificial viscosity for the temperature or
       * composition on each local cell.
       * @param viscosity_per_cell Output vector
       * @param advection_field Determines whether this variable should select
       * the temperature field or a compositional field.
       * @param skip_interior_cells A boolean flag. If set to true the function
       * will only compute the artificial viscosity in cells at boundaries.
       */
      template <typename T>
      void get_artificial_viscosity (Vector<T> &viscosity_per_cell,
                                     const AdvectionField &advection_field,
                                     const bool skip_interior_cells = false) const;

      /**
       * Adjust the pressure variable (which is only determined up to
       * a constant by the equations, though its value may enter
       * traction boundary conditions) by adding a constant to it in
       * such a way that the pressure on the surface or within the
       * entire volume has a known average value. The point of this
       * function is that the pressure that results from solving the
       * linear system may not coincide with what we think of as the
       * "physical pressure"; in particular, we typically think of the
       * pressure as zero (on average) along the surface because it is
       * the sum of the hydrostatic and dynamic pressure, where the
       * former is thought of as zero along the surface. This function
       * therefore converts from the "mathematical" pressure to the
       * "physical" pressure so that all following postprocessing
       * steps can use the latter.
       *
       * Whether the pressure should be normalized based on the
       * surface or volume average is decided by a parameter in the
       * input file.
       *
       * @note This function is called after setting the initial
       * pressure field in compute_initial_pressure_field() and at the end
       * of solve_stokes(). This makes sense because these are exactly the
       * places where the pressure is modified or re-computed.
       *
       * @return This function returns the pressure adjustment by value.
       * This is so that its negative can later be used again in
       * denormalize_pressure().
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      double normalize_pressure(LinearAlgebra::BlockVector &vector) const;

      /**
       * Invert the action of the normalize_pressure() function above. This
       * means that we move from a pressure that satisfies the pressure
       * normalization (e.g., has a zero average pressure, or a zero average
       * surface pressure) to one that does not actually satisfy this
       * normalization, and this doesn't seem to make sense because we are
       * not interested in such a pressure.
       *
       * Indeed, this function is only called at the very beginning of
       * solve_stokes() before we compute the initial (linear) residual
       * of the linear system $Ax=b$ that corresponds to the Stokes system,
       * where $x$ is the variable for which the pressure is adjusted
       * back to the "wrong" form. Because stokes_system() calls
       * normalize_pressure() at the end of its operations, no such
       * "wrong" pressure ever escapes the realm of solve_stokes(). The
       * "wrong" pressure is then used for two purposes in that function:
       * (i) To compute the initial Stokes residual, which makes sense
       * because it can only be zero (if we solved the same linear system
       * twice) if we re-use the exact same pressure as we got from the
       * previous solve -- i.e., before we called normalize_pressure()
       * at the end of the solve. (ii) To initialize the solution vector
       * before calling the GMRES solver, which also makes sense because
       * the best guess vector for GMRES is the one that had previously
       * come out of GMRES, namely the one on which we later called
       * normalize_pressure().
       *
       * This function modifies @p vector in-place.
       *
       * @note The adjustment made in this function is done using the
       * negative of the @p pressure_adjustment function argument that
       * would typically have been computed and returned by the
       * normalize_pressure() function. This value is typically stored in
       * the member variable @p last_pressure_normalization_adjustment,
       * but the current function doesn't read this variable but instead
       * gets the adjustment variable from the given argument.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void denormalize_pressure(const double                      pressure_adjustment,
                                LinearAlgebra::BlockVector       &vector) const;

      /**
       * Apply the bound preserving limiter to the discontinuous Galerkin solutions:
       * i.e., given two fixed upper and lower bound [min, max], after applying the limiter,
       * the discontinuous Galerkin solution will stay in the prescribed bounds.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void apply_limiter_to_dg_solutions (const AdvectionField &advection_field);

      /**
       * Compute the unique support points for the advection fields @p advection_fields.
       * The support points are collected by taking the union of the support points
       * of all given advection fields and filtering out duplicate points. The resulting
       * set of points is written into @p unique_support_points. @p support_point_index_by_field
       * is a vector of vectors that contains the indices of the support points for each field.
       * I.e. support_point_index_by_field[i][j] contains the j-th support point
       * for the i-th field in @p advection_fields and
       * unique_support_points[support_point_index_by_field[i][j]] is its location.
       *
       * For the common case that all fields have the same support points, each vector in
       * @p support_point_index_by_field will contain the same indices for all fields.
       *
       * Note that existing content of @p unique_support_points and @p support_point_index_by_field
       * will be overwritten in this function.
       */
      void compute_unique_advection_support_points (const std::vector<AdvectionField> &advection_fields,
                                                    std::vector<Point<dim>> &unique_support_points,
                                                    std::vector<std::vector<unsigned int>> &support_point_index_by_field) const;

      /**
       * Compute the reactions in case of operator splitting:
       * Using the current solution vector, this function makes a number of time
       * steps determined by the size of the reaction time step, and solves a
       * system of coupled ordinary differential equations for the reactions between
       * compositional fields and temperature in each of them. To do that, is uses
       * the reaction rates outputs from the material and heating models used in
       * the computation. The solution vector is then updated with the new values
       * of temperature and composition after the reactions.
       *
       * As the ordinary differential equation in any given point is independent
       * from the solution at all other points, we do not have to assemble a matrix,
       * but just need to loop over all node locations for the temperature and
       * compositional fields and compute the update to the solution.
       *
       * The function also updates the old solution vectors with the reaction update
       * so that the advection time stepping scheme will have the correct field terms
       * for the right-hand side when assembling the advection system.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void compute_reactions ();

      /**
       * Update the indicated block of the solution vector with the
       * corresponding block of the handed over @p distributed_vector. Also
       * update reaction_vector with the corresponding block of @p
       * distributed_reaction_vector.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void update_solution_vectors_with_reaction_results (const unsigned int block_index,
                                                          const LinearAlgebra::BlockVector &distributed_vector,
                                                          const LinearAlgebra::BlockVector &distributed_reaction_vector);

      /**
       * Initialize the current linearization point vector from the old
       * solution vector(s). Depending on the time of the call this
       * can be simply a copy of the solution of the last timestep,
       * or a linear extrapolation of old and old_old timestep to
       * the new timestep.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void initialize_current_linearization_point ();

      /**
       * Interpolate material model outputs onto an advection field (temperature
       * or composition). For the field identified by the AdvectionField @p adv_field, this function
       * asks the material model to fill a MaterialModel::MaterialModelOutputs
       * object that has an attached MaterialModel::PrescribedFieldOutputs
       * "additional outputs" object (for composition) or an attached
       * MaterialModel::PrescribedTemperatureOutputs (for temperature).
       * The MaterialModel::MaterialModelInputs
       * object passed to the material model then contains the support points of
       * the advection field, thereby allowing the outputs to be
       * interpolated to the finite element space and consequently into the
       * solution vector.
       * This is useful for advection fields whose advection mode is set
       * to Parameters::AdvectionFieldMethod::prescribed_field or
       * Parameters::AdvectionFieldMethod::prescribed_field_with_diffusion.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void interpolate_material_output_into_advection_field (const std::vector<AdvectionField> &adv_field);


      /**
       * Interpolate the given function onto the velocity FE space and write
       * it into the given vector.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void interpolate_onto_velocity_system(const TensorFunction<1,dim> &func,
                                            LinearAlgebra::Vector &vec) const;


      /**
       * Add constraints to the given @p constraints object that are required
       * for unique solvability of the velocity block based on the nullspace
       * removal settings.
       *
       * This method will add a zero Dirichlet constraint for the first
       * velocity unknown in the domain for each velocity component, which is
       * later being processed for translational or linear momentum removal.
       * This avoids breakdowns of the linear solvers that otherwise occurred
       * in some instances.
       *
       * @note: Rotational modes are currently not handled and don't appear to
       * require constraints so far.
       */
      void setup_nullspace_constraints(AffineConstraints<double> &constraints);


      /**
       * Eliminate the nullspace of the velocity in the given vector. Both
       * vectors are expected to contain the current solution.
       *
       * @param solution The locally relevant vector for the whole
       * finite element, this vector will be filled at the end.
       * @param distributed_stokes_solution only contains velocity and pressure and
       * only locally owned elements.
       *
       * This function is implemented in
       * <code>source/simulator/nullspace.cc</code>.
       */
      void remove_nullspace(LinearAlgebra::BlockVector &solution,
                            LinearAlgebra::BlockVector &distributed_stokes_solution) const;

      /**
       * Compute the angular momentum and other rotation properties
       * of the velocities in the given solution vector.
       *
       * @param use_constant_density determines whether to use a constant
       * density (which corresponds to computing a net rotation instead of net
       * angular momentum).
       * @param solution Solution vector to compute the properties for.
       * @param limit_to_top_faces allows to only compute the net angular momentum
       * (or net rotation) of the top surface.
       *
       * This function is implemented in
       * <code>source/simulator/nullspace.cc</code>.
       */
      RotationProperties<dim>
      compute_net_angular_momentum(const bool use_constant_density,
                                   const LinearAlgebra::BlockVector &solution,
                                   const bool limit_to_top_faces = false) const;

      /**
       * Remove the angular momentum of the given vector
       *
       * @param use_constant_density determines whether to use a constant
       * density (which corresponds to removing a net rotation instead of net
       * angular momentum).
       * @param relevant_dst locally relevant vector for the whole FE, will be
       * filled at the end.
       * @param tmp_distributed_stokes only contains velocity and pressure.
       * @param limit_to_top_faces allows to only remove the net angular momentum
       * (or net rotation) of the top surface. This can be useful to compare surface
       * motions against plate reconstructions in no net rotation reference frames.
       *
       * This function is implemented in
       * <code>source/simulator/nullspace.cc</code>.
       */
      void remove_net_angular_momentum(const bool use_constant_density,
                                       LinearAlgebra::BlockVector &relevant_dst,
                                       LinearAlgebra::BlockVector &tmp_distributed_stokes,
                                       const bool limit_to_top_faces = false) const;

      /**
       * Offset the boundary id of all faces located on an outflow boundary
       * by a fixed value given by the input parameter @p boundary_id_offset.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void replace_outflow_boundary_ids(const unsigned int boundary_id_offset,
                                        const bool is_composition,
                                        const unsigned int composition_index);

      /**
       * Undo the offset of the boundary ids done in replace_outflow_boundary_ids
       * by resetting all boundary ids to their original value.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void restore_outflow_boundary_ids(const unsigned int boundary_id_offset);

      /**
       * Remove the linear momentum of the given vector
       *
       * @param use_constant_density determines whether to use a constant
       * density (which corresponds to removing a net translation instead of
       * net linear momentum).
       * @param relevant_dst locally relevant vector for the whole FE, will be
       * filled at the end.
       * @param tmp_distributed_stokes only contains velocity and pressure.
       *
       * This function is implemented in
       * <code>source/simulator/nullspace.cc</code>.
       */
      void remove_net_linear_momentum(const bool use_constant_density,
                                      LinearAlgebra::BlockVector &relevant_dst,
                                      LinearAlgebra::BlockVector &tmp_distributed_stokes) const;

      /**
       * Compute the maximal velocity throughout the domain. This is needed to
       * compute the size of the time step.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      double get_maximal_velocity (const LinearAlgebra::BlockVector &solution) const;

      /**
       * Compute the variation (i.e., the difference between maximal and
       * minimal value) of the entropy $(T-\bar T)^2$ where $\bar T$ is the
       * average temperature throughout the domain given as argument to this
       * function.
       *
       * This function is used in computing the artificial diffusion
       * stabilization term.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      double get_entropy_variation (const double average_field,
                                    const AdvectionField &advection_field) const;

      /**
       * Compute the minimal and maximal temperature throughout the domain from
       * a solution vector extrapolated from the previous time steps. This is
       * needed to compute the artificial diffusion stabilization terms.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      std::pair<double,double>
      get_extrapolated_advection_field_range (const AdvectionField &advection_field) const;

      /**
       * Exchange coarsen/refinement flags set between processors so that
       * we have the correct settings on all ghost cells.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       *
       */
      void exchange_refinement_flags();


      /**
       * Check if timing output should be written in this timestep, and if so
       * write it.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       *
       */
      void maybe_write_timing_output () const;

      /**
       * Check if a checkpoint should be written in this timestep. If so create
       * one. Returns whether a checkpoint was written.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      bool maybe_write_checkpoint (const std::time_t last_checkpoint_time,
                                   const bool        force_writing_checkpoint);

      /**
       * Check if we should do an initial refinement cycle in this timestep.
       * This will only be checked in timestep 0, afterwards the variable
       * pre_refinement_step variable is invalidated, and this function will
       * return without doing refinement.
       * An initial refinement cycle is different from a regular one,
       * because time is not increased. Instead the same timestep is solved
       * using the new mesh.
       * Therefore, only output timing information and postprocessor output
       * if required in the input file. But always output statistics (to have
       * a history of the number of cells in the statistics file).
       * This function returns whether an initial refinement was done.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      bool maybe_do_initial_refinement (const unsigned int max_refinement_level);

      /**
       * Check if refinement is requested in this timestep. If so: Refine mesh.
       * The @p max_refinement_level might be increased from this time on
       * if this is an additional refinement cycle.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void maybe_refine_mesh (const double new_time_step,
                              unsigned int &max_refinement_level);

      /**
       * Advance the current time by the given @p step_size and update the
       * solution vectors as needed.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void advance_time (const double step_size);

      /**
       * Compute the artificial diffusion coefficient value on a cell given
       * the values and gradients of the solution passed as arguments.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      double
      compute_viscosity(internal::Assembly::Scratch::AdvectionSystem<dim> &scratch,
                        const double                        global_u_infty,
                        const double                        global_field_variation,
                        const double                        average_field,
                        const double                        global_entropy_variation,
                        const double                        cell_diameter,
                        const AdvectionField     &advection_field) const;

      /**
       * Compute the residual of one advection equation to be used for the
       * artificial diffusion coefficient value on a cell given the values and
       * gradients of the solution passed as arguments.
       *
       * This function is implemented in
       * <code>source/simulator/assembly.cc</code>.
       */
      void
      compute_advection_system_residual(internal::Assembly::Scratch::AdvectionSystem<dim> &scratch,
                                        const double                        average_field,
                                        const AdvectionField               &advection_field,
                                        double                             &max_residual,
                                        double                             &max_velocity,
                                        double                             &max_density,
                                        double                             &max_specific_heat,
                                        double                             &conductivity) const;

      /**
       * Return whether the Stokes matrix depends on the values of the
       * solution at the previous time step. This is the case is the
       * coefficients that appear in the matrix (i.e., the viscosity and, in
       * the case of a compressible model, the density) depend on the
       * solution.
       *
       * This function exists to ensure that the Stokes matrix is rebuilt in
       * time steps where it may have changed, while we want to save the
       * effort of rebuilding it whenever we don't need to.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      bool
      stokes_matrix_depends_on_solution () const;

      /**
       * Return whether to the best of our knowledge the A block of the
       * Stokes system is symmetric. This is the case for most models, except
       * if additional non-symmetric terms are added by special assemblers
       * (e.g., the free surface stabilization term).
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      bool
      stokes_A_block_is_symmetric () const;

      /**
       * This function checks that the user-selected formulations of the
       * equations are consistent with the other inputs. If an incorrect
       * selection is detected it throws an exception. It for example assures that
       * correct heating terms are selected, and the material model supports
       * the selection of the mass conservation formulation (e.g. incompressible)).
       * If the parameter 'parameters.formulation' is set to 'custom'
       * it only ensures very basic consistency.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void
      check_consistency_of_formulation ();

      /**
       * This function checks if the default solver and/or material
       * averaging were selected and if so, determines the appropriate
       * solver and/or averaging option.
       */
      void
      select_default_solver_and_averaging ();

      /**
       * This function checks that the user-selected boundary conditions do not
       * contain contradictions. If an incorrect selection is detected it
       * throws an exception. This for example assures that not both velocity
       * and traction boundary conditions are prescribed at the same boundary,
       * and that no boundary temperatures are prescribed at a periodic boundary.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void
      check_consistency_of_boundary_conditions () const;

      /**
       * Computes the initial Newton residual.
       */
      double
      compute_initial_newton_residual ();

      /**
       * This function computes the Eisenstat Walker linear tolerance used for the Newton iterations
       * in the `iterated Advection and Newton Stokes' and `single Advection, iterated Newton Stokes' solver schemes.
       * The Eisenstat and Walker (1996) method is used for determining the linear tolerance of
       * the iteration after the first iteration. The paper gives two preferred choices of computing
       * this tolerance. Both choices are implemented here with the suggested parameter values and
       * safeguards.
       */
      double
      compute_Eisenstat_Walker_linear_tolerance(const bool EisenstatWalkerChoiceOne,
                                                const double maximum_linear_stokes_solver_tolerance,
                                                const double linear_stokes_solver_tolerance,
                                                const double stokes_residual,
                                                const double newton_residual,
                                                const double newton_residual_old);

      /**
       * This function is called at the end of each time step and writes the
       * statistics object that contains data like the current time, the
       * number of linear solver iterations, and whatever the postprocessors
       * have generated, to disk.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      void output_statistics();

      /**
       * This routine computes the initial (nonlinear) Stokes residual that is
       * needed as a convergence criterion in models with solver schemes that do
       * nonlinear iterations. We calculate it in the same way as the tolerance for the linear
       * solver, using the norm of the pressure RHS for the pressure part and a
       * residual with zero velocity for the velocity part to get the part of
       * the RHS not balanced by the static pressure.
       *
       * This function is implemented in
       * <code>source/simulator/helper_functions.cc</code>.
       */
      double
      compute_initial_stokes_residual();

      /**
       * @}
       */

      /**
       * @name Variables that have to do with input, output, parallel
       * communication and interfacing with other parts of the program.
       * @{
       */
      Parameters<dim>                     parameters;

      /**
       * Unique pointer for an instance of the MeltHandler. This way,
       * if we do not need the machinery for doing melt stuff, we do
       * not even allocate it.
       */
      std::unique_ptr<MeltHandler<dim>> melt_handler;

      /**
       * Unique pointer for an instance of the NewtonHandler. This way,
       * if we do not need the machinery for doing Newton stuff, we do
       * not even allocate it.
       */
      std::unique_ptr<NewtonHandler<dim>> newton_handler;

      SimulatorSignals<dim>               signals;

      const IntermediaryConstructorAction post_signal_creation;

      /**
       * Unique pointer for an instance of the VolumeOfFluidHandler. This way,
       * if we do not need the machinery for doing volume_of_fluid stuff, we do
       * not even allocate it.
       *
       * Located here due to needing signals access
       */
      std::unique_ptr<VolumeOfFluidHandler<dim>> volume_of_fluid_handler;

      Introspection<dim>                  introspection;


      MPI_Comm                            mpi_communicator;

      /**
       * This stream will log into the file output/log.txt (used automatically
       * by pcout).
       */
      std::ofstream log_file_stream;

      using TeeDevice = boost::iostreams::tee_device<std::ostream, std::ofstream>;
      using TeeStream = boost::iostreams::stream<TeeDevice>;

      TeeDevice iostream_tee_device;
      TeeStream iostream_tee_stream;

      /**
       * Output stream for logging information. Will only output on processor
       * 0.
       */
      ConditionalOStream                  pcout;

      /**
       * An object that stores a bunch of statistics such as the number of
       * linear solver iterations, the time corresponding to each time step,
       * etc, as well as whatever the various postprocessors want to put into
       * it.
       *
       * This variable is written to disk after every time step, by the
       * Simulator::output_statistics() function.
       */
      TableHandler                        statistics;

      /**
       * The following two variables keep track which parts of the statistics
       * object have already been written. This is because the TableHandler
       * class has no way to keep track what it has already written, and so
       * we can not just append the last row of the table to the output
       * file. Rather, we keep track how many bytes we already wrote,
       * and a hash of what they contained, and if these so-many bytes have
       * not changed between the previous and current write operation, then
       * we only open the file in 'append' mode to add the new bytes from the
       * last row. If what we would write now has changed from what we wrote
       * back then in the first so-many bytes (e.g., because column widths of
       * the table have changed), then we just replace the previous file by
       * the current table contents in their entirety.
       */
      std::size_t                         statistics_last_write_size;
      std::size_t                         statistics_last_hash;

      mutable TimerOutput                 computing_timer;

      /**
       * A timer used to track the current wall time since the
       * last snapshot (or since the program started).
       */
      Timer wall_timer;

      /**
       * The total wall time that has elapsed up to the last snapshot
       * that was created.
       */
      double total_walltime_until_last_snapshot;

      /**
       * Checkpointing happens in rotating folders /restart/01/, /restart/02/,
       * etc.. This variable holds the last index used and as such should
       * contain the last valid checkpoint written.
       */
      unsigned int last_checkpoint_id;

      /**
       * In output_statistics(), where we output the statistics object above,
       * we do the actual writing on a separate thread. This variable is the
       * handle we get for this thread so that we can wait for it to finish,
       * either if we want to write the statistics object for the next thread,
       * or if we want to terminate altogether.
       */
      std::thread                         output_statistics_thread;

      /**
       * @}
       */

      /**
       * @name Variables that describe the physical setup of the problem
       * @{
       */
      const std::shared_ptr<InitialTopographyModel::Interface<dim>>          initial_topography_model;
      const std::unique_ptr<GeometryModel::Interface<dim>>                   geometry_model;
      const IntermediaryConstructorAction                                    post_geometry_model_creation_action;
      const std::unique_ptr<MaterialModel::Interface<dim>>                   material_model;
      const std::unique_ptr<GravityModel::Interface<dim>>                    gravity_model;

      BoundaryTemperature::Manager<dim>                                      boundary_temperature_manager;
      BoundaryConvectiveHeating::Manager<dim>                                boundary_convective_heating_manager;
      BoundaryComposition::Manager<dim>                                      boundary_composition_manager;
      const std::unique_ptr<PrescribedStokesSolution::Interface<dim>>        prescribed_stokes_solution;

      /**
       * The following two variables are pointers to objects that describe
       * the initial temperature and composition values. The Simulator
       * class itself releases these pointers once they are no longer
       * needed, somewhere during the first time step once it is known
       * that they are no longer needed. However, plugins can have their
       * own shared pointers to these objects, and the lifetime of the
       * objects pointed to is then until the last of these plugins
       * gets deleted.
       */
      std::shared_ptr<InitialTemperature::Manager<dim>>                      initial_temperature_manager;
      std::shared_ptr<InitialComposition::Manager<dim>>                      initial_composition_manager;

      const std::unique_ptr<AdiabaticConditions::Interface<dim>>             adiabatic_conditions;
#ifdef ASPECT_WITH_WORLD_BUILDER
      /**
       * A pointer to the WorldBuilder object. Like the
       * `initial_temperature_manager` and
       * `initial_composition_manager` objects above, the Simulator
       * object itself releases this pointer at the end of the
       * initialization process (right after releasing the
       * two mentioned initial condition objects). If a part of
       * the plugin system still needs the world builder object
       * after this point, it needs to keep its own shared pointer
       * to it.
       */
      std::shared_ptr<WorldBuilder::World>                      world_builder;
#endif
      BoundaryVelocity::Manager<dim>                            boundary_velocity_manager;
      BoundaryTraction::Manager<dim>                            boundary_traction_manager;
      const std::unique_ptr<BoundaryHeatFlux::Interface<dim>>   boundary_heat_flux;

      /**
       * @}
       */
      /**
       * @name Variables that describe the time discretization
       * @{
       */
      double                                                    time;
      double                                                    time_step;
      double                                                    old_time_step;
      unsigned int                                              timestep_number;
      unsigned int                                              pre_refinement_step;
      unsigned int                                              nonlinear_iteration;
      unsigned int                                              nonlinear_solver_failures;
      /**
       * @}
       */

      /**
       * @name Variables related to simulation time stepping
       * @{
       */
      TimeStepping::Manager<dim>                                time_stepping_manager;
      /**
       * @}
       */

      /**
       * @name Variables for doing lateral averaging
       * @{
       */
      LateralAveraging<dim>                                     lateral_averaging;
      /**
       * @}
       */

      /**
       * @name Variables that describe the spatial discretization
       * @{
       */
      parallel::distributed::Triangulation<dim>                 triangulation;
      double                                                    global_Omega_diameter;
      double                                                    global_volume;

      MeshRefinement::Manager<dim>                              mesh_refinement_manager;
      HeatingModel::Manager<dim>                                heating_model_manager;

      /**
       * Pointer to the Mapping object used by the finite elements when
       * going from the reference cell to the cell in the computational
       * domain. We use a pointer since different mapping objects may
       * be useful. In particular, when the mesh is deformable we use
       * a MappingQ1Eulerian object to describe the mesh deformation,
       * swapping it in for the original MappingQ or MappingCartesian object.
       */
      std::unique_ptr<Mapping<dim>>                             mapping;

      const FESystem<dim>                                       finite_element;

      DoFHandler<dim>                                           dof_handler;

      Postprocess::Manager<dim>                                 postprocess_manager;

      /**
       * The managers holding different sets of particles
       */
      std::vector<Particle::Manager<dim>>                       particle_managers;

      /**
       * Constraint objects. The first of these describes all constraints that
       * are not time dependent (e.g., hanging nodes, no-normal-flux
       * constraints), whereas the second one is initialized at the top of
       * every time step by copying from the first and then adding to it
       * constraints that are time dependent (e.g., time dependent velocity or
       * temperature boundary conditions).
       *
       * 'constraints' is computed in setup_dofs(), 'current_constraints' is
       * done in compute_current_constraints().
       */
      AffineConstraints<double>                                 constraints;
      AffineConstraints<double>                                 current_constraints;

      /**
       * A place to store the latest correction computed by normalize_pressure().
       * We store this so we can undo the correction in denormalize_pressure().
       */
      double                                                    last_pressure_normalization_adjustment;

      /**
       * Scaling factor for the pressure as explained in the
       * Kronbichler/Heister/Bangerth paper to ensure that the linear system
       * that results from the Stokes equations is well conditioned.
       */
      double                                                    pressure_scaling;

      /**
       * A variable that determines whether we need to do the correction of
       * the Stokes right hand side vector to ensure that the average
       * divergence is zero. This is necessary for compressible models, but
       * only if there are no in/outflow boundaries.
       */
      bool                           do_pressure_rhs_compatibility_modification;

      /**
       * @}
       */


      /**
       * @name Variables that describe the linear systems and solution vectors
       * @{
       */

      /**
       * An object that contains the entries of the system matrix. It
       * has a size equal to the total number of degrees of freedom,
       * but since we typically do not solve for all variables at
       * once, the content of the matrix at any given time is only
       * appropriate for the part of the system we are currently
       * solving.
       */
      LinearAlgebra::BlockSparseMatrix                          system_matrix;

      /**
       * This vector is used for the weighted BFBT preconditioner. It
       * stores the inverted lumped velocity mass matrix.
       */
      LinearAlgebra::BlockVector                                inverse_lumped_mass_matrix;

      /**
       * An object that contains the entries of preconditioner
       * matrices for the system matrix. It has a size equal to the
       * total number of degrees of freedom, but is only used for the
       * Stokes system (that's the only part of the system where we
       * use a matrix for preconditioning that is different from the
       * matrix we solve). Consequently, the blocks in rows and
       * columns corresponding to temperature or compositional fields
       * are left empty when building the sparsity pattern of this
       * matrix in the Simulator::setup_system_preconditioner()
       * function.
       */
      LinearAlgebra::BlockSparseMatrix                          system_preconditioner_matrix;

      LinearAlgebra::BlockVector                                solution;
      LinearAlgebra::BlockVector                                old_solution;
      LinearAlgebra::BlockVector                                old_old_solution;
      LinearAlgebra::BlockVector                                system_rhs;

      LinearAlgebra::BlockVector                                current_linearization_point;

      // only used if is_compressible()
      LinearAlgebra::BlockVector                                pressure_shape_function_integrals;

      // only used if operator split is enabled
      LinearAlgebra::BlockVector                                operator_split_reaction_vector;



      std::unique_ptr<LinearAlgebra::PreconditionAMG>           Amg_preconditioner;
      std::unique_ptr<LinearAlgebra::PreconditionBase>          Mp_preconditioner;

      bool                                                      rebuild_sparsity_and_matrices;
      bool                                                      rebuild_stokes_matrix;
      bool                                                      assemble_newton_stokes_matrix;
      bool                                                      assemble_newton_stokes_system;
      bool                                                      rebuild_stokes_preconditioner;

      /**
       * @}
       */

    private:

      /**
       * Unique pointer for an instance of the MeshDeformationHandler. this way,
       * if we do not need the machinery for doing mesh deformation stuff, we do
       * not even allocate it.
       */
      std::unique_ptr<MeshDeformation::MeshDeformationHandler<dim>> mesh_deformation;

      /**
       * Unique pointer for the matrix-free Stokes solver
       */
      std::unique_ptr<StokesMatrixFreeHandler<dim>> stokes_matrix_free;

      /**
       * Unique pointer for the direct Stokes solver
       */
      std::unique_ptr<StokesSolver::Direct<dim>> stokes_direct;


      friend class boost::serialization::access;
      friend class SimulatorAccess<dim>;
      friend class MeshDeformation::MeshDeformationHandler<dim>;
      friend class VolumeOfFluidHandler<dim>;
      friend class StokesMatrixFreeHandler<dim>;
      template <int dimension, int velocity_degree> friend class StokesMatrixFreeHandlerLocalSmoothingImplementation;
      friend struct Parameters<dim>;
  };
}


#endif
