/*
  Copyright (C) 2019 - 2020 by the authors of the ASPECT code.

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


#ifndef _aspect_time_stepping_interface_h
#define _aspect_time_stepping_interface_h

#include <aspect/global.h>
#include <aspect/simulator_access.h>
#include <aspect/termination_criteria/interface.h>

namespace aspect
{
  using namespace dealii;

  /**
   * A namespace containing all class related to the time stepping plugin system.
   */
  namespace TimeStepping
  {

    /**
     * This enum describes the possible reactions by the time stepping plugins after the
     * computation of each time step.
     *
     * Note: the ordering of the options is crucial, as the Manager will use the minimum
     * of the values given by all active plugins as the reaction to take. This means the
     * first entry has highest priority.
     */
    enum class Reaction
    {
      /**
       * Initiate mesh refinement and go back in time to repeat the last timestep.
       */
      refine_and_repeat_step,
      /**
       * Go back in time to repeat the last timestep.
       */
      repeat_step,
      /**
       * Initiate mesh refinement and continue to the next timestep.
       */
      refine_and_advance,
      /**
       * Continue to the next timestep. The default action to take.
       */
      advance
    };

    /**
     * Information passed to Interface::determine_reaction with information
     * about the time step.
     */
    struct TimeStepInfo
    {
      /**
       * The proposed time step size for the next step as computed by querying
       * plugins and applying other logic.
       */
      double next_time_step_size;

      /**
       * If true, a termination criterion decided to shorten the time step
       * size.
       */
      bool reduced_by_termination_plugin;
    };

    /**
     * A base class for parameterizations of the time stepping models.
     *
     * @ingroup TimeStepping
     */
    template <int dim>
    class Interface
    {
      public:
        /**
         * Destructor. Made virtual to enforce that derived classes also have
         * virtual destructors.
         */
        virtual ~Interface() = default;

        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after
         * the SimulatorAccess (if applicable) is initialized.
         */
        virtual
        void
        initialize ();

        /**
         * A function that is called at the beginning of each time step. The
         * default implementation of the function does nothing, but derived
         * classes that need more elaborate setups for a given time step may
         * overload the function.
         */
        virtual
        void
        update ();

        /**
         * Execute the logic of the plugin.
         *
         * This is called after every time step to determine
         * a) What to do (advance, repeat, etc.), see the Reaction enum.
         * b) What timestep size to use.
         *
         */
        virtual
        double
        execute() = 0;

        /**
         * Determine what we want with the simulation to happen next: advance,
         * repeat, refinement, etc.. The second return value is the time step
         * size to take in case the plugin requests a repeated time step.
         *
         * The argument @p info contains information like the step size that
         * would be taken in this time step (determined as the minimum of the
         * return value of execute() from all plugins).
         *
         * The default implementation of this function will always advance
         * to the next time step.
         */
        virtual
        std::pair<Reaction, double>
        determine_reaction(const TimeStepInfo &info);

        /**
         * Declare the parameters this class takes through input files. The
         * default implementation of this function does not describe any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * The default implementation of this function does not read any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);
    };


    /**
     * A class to handle computation of the next time step (as desired by the user) and
     * checking if the simulation is finished.
     */
    template <int dim>
    class Manager : public SimulatorAccess<dim>
    {
      public:
        /**
         * Override initialize_simulator() so that we can also initialize the contained
         * termination_manager.
         */
        virtual void initialize_simulator (const Simulator<dim> &simulator_object) override;


        /**
         * Update the current state and determine what needs to happen based on the
         * last computed solution (see functions below). This computes the size of
         * the next time step potentially taking into account the current solution
         * (convection time step, conduction time step), settings from parameters,
         * and termination criteria (to hit the end time exactly).
         */
        void update();

        /**
         * Return the next step size as computed from update().
         */
        double get_next_time_step_size() const;

        /**
         * If true, a plugin requested to redo the last computed time step. Updated
         * when calling update().
         */
        bool should_repeat_time_step() const;

        /**
          * If true, execute a mesh refinement step now (potentially before repeating
          * the current time step).
          */
        bool should_refine_mesh() const;

        /**
         * If true, the simulator should perform a checkpoint before terminating.
         */
        bool need_checkpoint_on_terminate() const;

        /**
         * Check if the simulation is ready to terminate successfully.
         */
        bool should_simulation_terminate_now() const;

        /**
         * Declare the parameters of all known termination criteria plugins,
         * as well as of ones this class has itself.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * This determines which termination criteria objects will be created;
         * then let these objects read their parameters as well.
         */
        void
        parse_parameters (ParameterHandler &prm);

        /**
         * For the current plugin subsystem, write a connection graph of all of the
         * plugins we know about, in the format that the
         * programs dot and neato understand. This allows for a visualization of
         * how all of the plugins that ASPECT knows about are interconnected, and
         * connect to other parts of the ASPECT code.
         *
         * @param output_stream The stream to write the output to.
         */
        static
        void
        write_plugin_graph (std::ostream &output_stream);


        /**
         * A function that is used to register time stepping model objects in such
         * a way that the Manager can deal with all of them without having to
         * know them by name. This allows the files in which individual
         * plugins are implemented to register these plugins, rather than also
         * having to modify the Manager class by adding the new plugin class.
         *
         * @param name A string that identifies the model
         * @param description A text description of what this model does and that
         * will be listed in the documentation of the parameter file.
         * @param declare_parameters_function A pointer to a function that can be
         * used to declare the parameters that this model wants to read
         * from input files.
         * @param factory_function A pointer to a function that can create an
         * object of this model.
         *
         * @ingroup TimeStepping
         */
        static
        void
        register_time_stepping_model (const std::string &name,
                                      const std::string &description,
                                      void (*declare_parameters_function) (ParameterHandler &),
                                      Interface<dim> *(*factory_function) ());

      private:

        /**
         * The current Reaction computed by update().
         */
        Reaction current_reaction;

        /**
         * The next time step size computed by update().
         */
        double next_time_step_size;

        /**
         * The minimum time step size specified by the user (in seconds).
         */
        double minimum_time_step_size;

        /**
         * Whether to do a final checkpoint before termination. This is
         * specified in the parameters.
         */
        bool do_checkpoint_on_terminate;

        /**
         * The termination manager keeps track of the termination plugins and we use
         * it to determine the time_step size in the final time step.
         */
        TerminationCriteria::Manager<dim> termination_manager;

        /**
         * A list of active plugins to determine time step sizes.
         */
        std::list<std::unique_ptr<Interface<dim> > > active_plugins;
    };

    /**
     * Given a class name, a name, and a description for the parameter file, register it with the
     * aspect::TimeStepping::Manager class.
     *
     * @ingroup TimeStepping
     */
#define ASPECT_REGISTER_TIME_STEPPING_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_TIME_STEPPING_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::TimeStepping::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::TimeStepping::Manager<2>::register_time_stepping_model, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::TimeStepping::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::TimeStepping::Manager<3>::register_time_stepping_model, \
                                name, description); \
  }

  }
}

#endif
