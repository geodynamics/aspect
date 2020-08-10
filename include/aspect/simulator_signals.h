/*
  Copyright (C) 2015 - 2020 by the authors of the ASPECT code.

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


#ifndef _aspect_simulator_signals_h
#define _aspect_simulator_signals_h

#include <aspect/global.h>
#include <aspect/simulator_access.h>
#include <aspect/parameters.h>

#include <deal.II/base/parameter_handler.h>

#include <boost/signals2.hpp>

namespace aspect
{
  namespace Assemblers
  {
    template <int dim>
    class Manager;

    template <int dim>
    class Interface;
  }

  /**
   * A class that collects the definition of signals that can be triggered
   * at different points in a computation. A signal is in essence an event
   * that is triggered whenever the program passes a certain point in a
   * computation. Parties interested in any of these signals can attach
   * "slots" to a signal. A slot is, in essence, a function that is called
   * whenever the signal is triggered. Multiple slots (or none) can be
   * attached to the same signal. To be as general as possible, slots are
   * not actually just pointers to functions, but std::function objects
   * that have a certain signature. Consequently, they can have much more
   * complicated types than just function pointers, such as objects with
   * an <code>operator()</code> or lambda functions.
   *
   * The documentation of each of the signals below indicates when
   * exactly it is called.
   *
   * @ingroup Simulator
   */
  template <int dim>
  struct SimulatorSignals
  {
    /**
     * A signal that is called before the list of finite element variables is
     * used to construct the Introspection class.
     *
     * The functions (slots) that can attach to this signal need to
     * take one argument: A std::vector of VariableDeclaration<dim>
     * representing the collection of finite element variables,
     * that can be modified and will be used to construct the
     * final finite element system later.
     */
    boost::signals2::signal<void (std::vector<VariableDeclaration<dim> > &)>
    edit_finite_element_variables;

    /**
     * A signal that is called before setting up the initial conditions.
     *
     * The functions (slots) that can attach to this signal need to take one
     * argument: A SimulatorAccess object that describes the simulator.
     */
    boost::signals2::signal<void (typename parallel::distributed::Triangulation<dim> &)>  pre_set_initial_state;

    /**
     * A signal that is called after setting up the initial conditions.
     *
     * The functions (slots) that can attach to this signal need to take one
     * argument: A SimulatorAccess object that describes the simulator.
     */
    boost::signals2::signal<void (const SimulatorAccess<dim> &)>  post_set_initial_state;

    /**
     * A signal that is called at the end of setting up the
     * constraints for the current time step. This allows to add
     * more constraints on degrees of freedom, for example to fix
     * the velocity at certain points.
     *
     * The functions (slots) that can attach to this signal need to
     * take two arguments: A SimulatorAccess object that
     * describes the simulator to act on, and an (output)
     * argument that indicates the constraints to be computed.
     */
    boost::signals2::signal<void (const SimulatorAccess<dim> &,
                                  AffineConstraints<double> &)>  post_constraints_creation;

    /**
     * A signal that is called at the start of setup_dofs(). This allows for
     * editing of the parameters struct on the fly (such as changing boundary
     * conditions) to give ASPECT different behavior in mid-run than it
     * otherwise would have.
     *
     * The functions that connect to this signal must take two arguments, a
     * SimulatorAccess object that describes the simulator, and an object of
     * type aspect::Parameters<dim>, which is the current parameters object
     * that the simulator is working with.
     */
    boost::signals2::signal<void (const SimulatorAccess<dim> &,
                                  Parameters<dim> &parameters)>  edit_parameters_pre_setup_dofs;

    /**
     * A signal that is called before every mesh_refinement. This signal
     * allows for registering functions that store data that is related to
     * mesh cells and needs to be transferred with the cells during the
     * repartitioning.
     *
     * The functions that connect to this signal must take a reference to a
     * parallel::distributed::Triangulation object as argument. This argument
     * will point to the triangulation used by the Simulator class.
     */
    boost::signals2::signal<void (typename parallel::distributed::Triangulation<dim> &)>  pre_refinement_store_user_data;

    /**
     * A signal that is called after every mesh_refinement.  This signal
     * allows for registering functions that load data related to mesh cells
     * and that was transferred with the cells during the repartitioning.
     *
     * The functions that connect to this signal must take a reference to a
     * parallel::distributed::Triangulation object as argument. This argument
     * will point to the triangulation used by the Simulator class.
     */
    boost::signals2::signal<void (typename parallel::distributed::Triangulation<dim> &)>  post_refinement_load_user_data;

    /**
     * A signal that is called before the computation of tangential boundary
     * conditions for which normal vectors are needed, i.e. calls to the
     * compute_no_normal_flux_constraints function for both the
     * velocity variable in the main simulator and the mesh velocity
     * variable in models with mesh deformation.
     *
     * The functions that connect to this signal must take a reference
     * to a parallel::distributed::Triangulation object as
     * argument. This argument will point to the triangulation used by
     * the Simulator class.
     */
    boost::signals2::signal<void (typename parallel::distributed::Triangulation<dim> &)>  pre_compute_no_normal_flux_constraints;

    /**
     * A signal that is called after the computation of tangential boundary
     * conditions for which normal vectors are needed, i.e. calls to the
     * compute_no_normal_flux_constraints function for both the
     * velocity variable in the main simulator and the mesh velocity
     * variable in models with mesh deformation.
     *
     * The functions that connect to this signal must take a reference
     * to a parallel::distributed::Triangulation object as
     * argument. This argument will point to the triangulation used by
     * the Simulator class.
     */
    boost::signals2::signal<void (typename parallel::distributed::Triangulation<dim> &)>  post_compute_no_normal_flux_constraints;

    /**
     * A signal that is called before the creation of every checkpoint.  This
     * signal allows for registering functions that store data related to mesh
     * cells.
     *
     * The functions that connect to this signal must take a reference to a
     * parallel::distributed::Triangulation object as argument. This argument
     * will point to the triangulation used by the Simulator class.
     */
    boost::signals2::signal<void (typename parallel::distributed::Triangulation<dim> &)>  pre_checkpoint_store_user_data;

    /**
     * A signal that is called after resuming from a checkpoint.  This signal
     * allows for registering functions that load data related to mesh cells
     * that was previously stored in the checkpoint.  Note that before calling
     * Triangulation::notify_ready_to_unpack() the function needs to call
     * register_attach_data() with the appropriate arguments to restore the
     * state of the triangulation.
     *
     * The functions that connect to this signal must take a reference to a
     * parallel::distributed::Triangulation object as argument. This argument
     * will point to the triangulation used by the Simulator class.
     */
    boost::signals2::signal<void (typename parallel::distributed::Triangulation<dim> &)>  post_resume_load_user_data;

    /**
     * A signal that is called at the beginning of the program. It
     * gives user extensions the ability to declare additional
     * parameters via the provided argument. User extensions connected to
     * this signal will likely also want to connect to the
     * parse_additional_parameters signal.
     *
     * The first argument to functions that connect to this signal
     * denotes the dimension in which ASPECT will be run. Functions
     * connected to this signal can declare additional runtime
     * parameters in the second argument.
     */
    static boost::signals2::signal<void (const unsigned int aspect_dim,
                                         ParameterHandler &prm)>  declare_additional_parameters;

    /**
     * A signal that is called at the beginning of the program, after reading
     * the input file. It gives user extensions the ability to read additional
     * parameters from the provided argument. The first argument indicates an
     * object that represents all of the other parameters (that have already
     * been parsed at this point).
     *
     * User extensions connected to this signal will likely also want to
     * connect to the declare_additional_parameters signal.
     */
    static boost::signals2::signal<void (const Parameters<dim> &,
                                         ParameterHandler &)>  parse_additional_parameters;

    /**
     * A signal that is fired when the iterative Stokes solver is done.
     * Parameters are a reference to the SimulatorAccess, the number of
     * preconditioner inner solver iterations for the S and A block of the
     * system, and two information objects that contain information
     * about the success of the solve, the number of outer GMRES iterations
     * and the residual history for the cheap and expensive solver phase.
     */
    boost::signals2::signal<void (const SimulatorAccess<dim> &,
                                  const unsigned int number_S_iterations,
                                  const unsigned int number_A_iterations,
                                  const SolverControl &solver_control_cheap,
                                  const SolverControl &solver_control_expensive)> post_stokes_solver;

    /**
     * A signal that is fired when the iterative advection solver is done.
     * Parameters are a reference to the SimulatorAccess, a bool indicating
     * whether the temperature field or a compositional field was solved,
     * a composition index that describes which compositional field
     * was solved, and an information object that contains information
     * about the number of iterations and history of residuals.
     */
    boost::signals2::signal<void (const SimulatorAccess<dim> &,
                                  const bool solved_temperature_field,
                                  const unsigned int compositional_index,
                                  const SolverControl &solver_control)> post_advection_solver;

    /**
     * A signal that is fired when the nonlinear solver scheme is done.
     * The signal parameter is an object that contains information
     * about the final state (failure/success), number of
     * iterations and history of residuals of the nonlinear solver.
     * If there is no nonlinear solver (only a single solve), the
     * SolverControl object will report a successful state, a single iteration
     * and a remaining residual of zero.
     */
    boost::signals2::signal<void (const SolverControl &)> post_nonlinear_solver;

    /**
     * A signal that is fired at the end of the set_assemblers() function that
     * allows modification of the assembly objects active in this simulation.
     */
    boost::signals2::signal<void (const SimulatorAccess<dim> &,
                                  aspect::Assemblers::Manager<dim> &)>
    set_assemblers;
  };


  namespace internals
  {
    /**
     * A namespace for some internal functions that have to do with how plugins
     * can register their slots with signals.
     */
    namespace SimulatorSignals
    {
      /**
       * Two functions that (in 2d and 3d) put a user-provided function onto a list
       * of functions that the Simulator object will later go through when
       * letting plugins connect their slots to signals.
       */
      void register_connector_function_2d (const std::function<void (aspect::SimulatorSignals<2> &)> &connector);
      void register_connector_function_3d (const std::function<void (aspect::SimulatorSignals<3> &)> &connector);

      /**
       * A function that is called by the Simulator object and that goes
       * through the list (with the corresponding dimension) created by the
       * previous pair of functions and call each of the user-provided
       * connector functions to let them register their slots with the
       * corresponding signals.
       */
      template <int dim>
      void call_connector_functions (aspect::SimulatorSignals<dim> &signals);
    }
  }


  /**
   * A macro that is used in user-provided plugins to register a function that
   * is called at the beginning of a simulation by a Simulator object. When called,
   * the provided function will receive a SimulatorSignals object that contains
   * signals to which one can subscribe.
   *
   * For technical reasons, the macro takes two arguments denoting functions for the
   * 2d and 3d cases. These can, for example, be the names of 2d and 3d
   * instantiations of the same template function.
   */
#define ASPECT_REGISTER_SIGNALS_CONNECTOR(connector_function_2d,connector_function_3d) \
  namespace ASPECT_REGISTER_SIGNALS_CONNECTOR \
  { \
    struct dummy_do_register \
    {          \
      dummy_do_register () \
      { \
        aspect::internals::SimulatorSignals::register_connector_function_2d (connector_function_2d); \
        aspect::internals::SimulatorSignals::register_connector_function_3d (connector_function_3d); \
      } \
    } dummy_variable; \
  }


  /**
   * A macro that is used to register a function that can be used to connect user
   * extension functions to the parameter-related signals declared in SimulatorSignals.
   *
   * In essence, this function simply registers a (global) function that is called
   * at the beginning of the program and that can be used to connect parameter
   * declaration and parsing functions to the signals listed above.
   */
#define ASPECT_REGISTER_SIGNALS_PARAMETER_CONNECTOR(connector_function) \
  namespace ASPECT_REGISTER_SIGNALS_PARAMETER_CONNECTOR_ ## connector_function \
  { \
    struct dummy_do_register_ ## connector_function \
    { \
      dummy_do_register_ ## connector_function () \
      {                \
        connector_function (); \
      }          \
    } dummy_variable_ ## classname; \
  }

}
#endif
