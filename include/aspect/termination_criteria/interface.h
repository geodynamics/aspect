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


#ifndef _aspect_termination_criteria_interface_h
#define _aspect_termination_criteria_interface_h

#include <aspect/global.h>
#include <aspect/plugins.h>

#include <memory>
#include <deal.II/base/parameter_handler.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  template <int dim> class Simulator;


  /**
   * A namespace for everything to do with determining criteria for
   * terminating the simulation gracefully. This does not include termination
   * due to errors.
   *
   * @ingroup TerminationCriteria
   */
  namespace TerminationCriteria
  {

    /**
     * This class declares the public interface of termination criteria
     * plugins. These plugins must implement a function that can be called
     * each time step to determine if specific criteria have been met that
     * mean the simulation should be ended gracefully.
     *
     * Access to the data of the simulator is granted by the @p protected
     * member functions of the SimulatorAccess class, i.e., classes
     * implementing this interface will in general want to derive from both
     * this Interface class as well as from the SimulatorAccess class if they
     * need to find out about the state of the simulation.
     *
     * @ingroup TerminationCriteria
     */
    template <int dim>
    class Interface : public Plugins::InterfaceBase
    {
      public:
        /**
         * Execute evaluation of the termination criterion.
         *
         * @return Whether to terminate the simulation (true) or continue
         * (false).
         *
         * @note For some kinds of termination plugins, it may be difficult to
         * ensure that all processors come to the same conclusion. An example
         * is one that terminates the simulation after a certain amount of CPU
         * time has expired, and this may be different from MPI process to MPI
         * process. To avoid difficult to deal with situations, the plugin
         * manager only requires that a single processor returns
         * <code>true</code> for a given plugin to record that the plugin
         * wants to terminate the simulation, even if the other processors
         * return false.
         */
        virtual
        bool
        execute () = 0;

        /**
         * Check for last time step and if so reduce the time step to user
         * specified end time
         *
         * @return Reduced or current time step size
         *
         * @note A reduced time step size may be returned for the last time
         * step. For all other time steps, the current time step size
         * (provided as argument) will be returned. The returned time step
         * size will be greater than zero, and less than or equal to the given
         * argument put into this function.
         */
        virtual double check_for_last_time_step (const double time_step) const;
    };






    /**
     * A class that manages all objects that provide functionality to specify
     * termination criteria.
     *
     * @ingroup TerminationCriteria
     */
    template <int dim>
    class Manager : public Plugins::ManagerBase<Interface<dim>>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Execute all of the termination criteria objects that have been
         * requested in the input file.
         *
         * The function returns a bool indicating whether the simulation should be
         * terminated.
         *
         * To avoid undefined situations, this function only requires that a
         * single processor's termination request comes back positive for a
         * given plugin. In other words, for each plugin selected, not all
         * processors need to return the same value: if only one of the says
         * that the simulation should be terminated, then this is enough.
         */
        virtual
        bool
        execute () const;

        /**
         * Check all of the termination criteria objects that have been
         * requested in the input file for criteria regarding last time step
         * and if so get the minimum of these values.
         *
         * @return Reduced or current time step size
         *
         * @note A reduced time step size may be returned for the last time
         * step. For all other time steps, the current time step size
         * (provided as argument) will be returned. The returned time step
         * size will be greater than zero, and less than or equal to the given
         * argument put into this function.
         */
        double check_for_last_time_step (const double time_step) const;

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
        parse_parameters (ParameterHandler &prm) override;

        /**
         * A function that is used to register termination criteria objects in
         * such a way that the Manager can deal with all of them without
         * having to know them by name. This allows the files in which
         * individual plugins are implement to register these plugins, rather
         * than also having to modify the Manager class by adding the new
         * termination criteria class.
         *
         * @param name The name under which this plugin is to be called in
         * parameter files.
         * @param description A text description of what this model does and
         * that will be listed in the documentation of the parameter file.
         * @param declare_parameters_function A pointer to a function that
         * declares the parameters for this plugin.
         * @param factory_function A pointer to a function that creates such a
         * termination criterion object and returns a pointer to it.
         */
        static
        void
        register_termination_criterion (const std::string &name,
                                        const std::string &description,
                                        void (*declare_parameters_function) (ParameterHandler &),
                                        std::unique_ptr<Interface<dim>> (*factory_function) ());

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
         * Exception.
         */
        DeclException1 (ExcTerminationCriteriaNameNotFound,
                        std::string,
                        << "Could not find entry <"
                        << arg1
                        << "> among the names of registered termination criteria objects.");
    };



    /**
     * Given a class name, a name, and a description for the parameter file
     * for a termination criterion object, register it with the
     * aspect::TerminationCriteria::Manager class.
     *
     * @ingroup TerminationCriteria
     */
#define ASPECT_REGISTER_TERMINATION_CRITERION(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_TERMINATION_CRITERION_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::TerminationCriteria::Interface<2>,classname<2>> \
    dummy_ ## classname ## _2d (&aspect::TerminationCriteria::Manager<2>::register_termination_criterion, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::TerminationCriteria::Interface<3>,classname<3>> \
    dummy_ ## classname ## _3d (&aspect::TerminationCriteria::Manager<3>::register_termination_criterion, \
                                name, description); \
  }
  }
}


#endif
