/*
  Copyright (C) 2015 by the authors of the ASPECT code.

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


#ifndef __aspect__simulator_signals_h
#define __aspect__simulator_signals_h

#include <aspect/global.h>
#include <aspect/simulator_access.h>
#include <aspect/parameters.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/lac/constraint_matrix.h>

#include <boost/signals2.hpp>

namespace aspect
{
  namespace internal
  {
    namespace Assembly
    {
      namespace Assemblers
      {
        template <int dim>
        class AssemblerBase;
      }

      template <int dim>
      struct AssemblerLists;
    }
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
   * an <code>operator()</code> or function calls treated with things
   * like std::bind.
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
                                  ConstraintMatrix &)>  post_constraints_creation;

    /**
    * A signal that is called at the start of setup_dofs().
    * This allows for editing of the parameters struct on the fly
    * (such as changing boundary conditions) to give Aspect different behavior
    * in mid-run than it otherwise would have.

    * The functions that connect to this signal must take two arguments,
    * a SimulatorAccess object that describes the simulator, and an object
    * of type aspect::Parameters<dim>, which is the current parameters
    * object that the simulator is working with.
    */
    boost::signals2::signal<void (const SimulatorAccess<dim> &,
                                  Parameters<dim> &parameters)>  edit_parameters_pre_setup_dofs;

    /**
    * A signal that is called before every mesh_refinement.
    * This signal allows for registering functions that store data
    * that is related to mesh cells and needs to be transferred with the cells
    * during the repartitioning.

    * The functions that connect to this signal must take a list of pairs of
    * a size and a callback function that has the signature of the
    * parallel::distributed::Triangulation::register_attach_data function.
    * The intention of the function should be to add another pair, describing
    * the amount of space per cell that needs to be send and a function
    * that is called for every cell by the Triangulation.
    */
    boost::signals2::signal<void (typename parallel::distributed::Triangulation<dim> &)>  pre_refinement_store_user_data;

    /**
    * A signal that is called after every mesh_refinement.
    * This signal allows for registering functions that load data
    * that is related to mesh cells and was transferred with the cells
    * during the repartitioning.

    * The functions that connect to this signal must take a list of callback
    * functions that have the signature of the
    * parallel::distributed::Triangulation::register_attach_data function.
    * The intention of the function should be to add another function,
    * that is called for every cell by the Triangulation.
    */
    boost::signals2::signal<void (typename parallel::distributed::Triangulation<dim> &)>  post_refinement_load_user_data;

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
      * A signal that is fired when the iterative Stokes solver is
      * done. Parameters are a reference to the SimulatorAccess, a bool
      * indicating success or failure, and a vector with linear residuals in
      * each solver step. If the solver switches from the cheap to the expensive
      * solver, a -1.0 is inserted.
      */
    boost::signals2::signal<void (const SimulatorAccess<dim> &,
                                  const bool,
                                  const std::vector<double> &)> post_stokes_solver;


    /**
      * A signal that is fired at the end of the set_assemblers() function that allows
      * modification of the assembly objects active in this simulation.
      */
    boost::signals2::signal<void (const SimulatorAccess<dim> &,
                                  aspect::internal::Assembly::AssemblerLists<dim> &,
                                  std::vector<std_cxx11::shared_ptr<internal::Assembly::Assemblers::AssemblerBase<dim> > > &)>
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
      void register_connector_function_2d (const std_cxx11::function<void (aspect::SimulatorSignals<2> &)> &connector);
      void register_connector_function_3d (const std_cxx11::function<void (aspect::SimulatorSignals<3> &)> &connector);

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
    int dummy_do_register () \
    { \
      aspect::internals::SimulatorSignals::register_connector_function_2d (connector_function_2d); \
      aspect::internals::SimulatorSignals::register_connector_function_3d (connector_function_3d); \
      return /* anything will do = */42; \
    } \
    \
    const int dummy_variable = dummy_do_register (); \
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
    int dummy_do_register_ ## connector_function () \
    { \
      connector_function (); \
      return /* anything will do = */42; \
    } \
    \
    const int dummy_variable_ ## classname = dummy_do_register_ ## connector_function (); \
  }

}
#endif
