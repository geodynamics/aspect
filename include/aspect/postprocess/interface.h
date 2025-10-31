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


#ifndef _aspect_postprocess_interface_h
#define _aspect_postprocess_interface_h

#include <aspect/global.h>
#include <aspect/plugins.h>
#include <aspect/simulator_access.h>

#include <memory>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/parameter_handler.h>

#include <boost/serialization/split_member.hpp>
#include <boost/core/demangle.hpp>

#include <typeinfo>


namespace aspect
{
  template <int dim> class Simulator;
  template <int dim> class SimulatorAccess;


  /**
   * A namespace for everything to do with postprocessing solutions every time
   * step or every few time steps.
   *
   * @ingroup Postprocessing
   */
  namespace Postprocess
  {

    /**
     * This class declares the public interface of postprocessors.
     * Postprocessors must implement a function that can be called at the end
     * of each time step to evaluate the current solution, as well as
     * functions that save the state of the object and restore it (for
     * checkpoint/restart capabilities).
     *
     * Access to the data of the simulator is granted by the @p protected
     * member functions of the SimulatorAccess class, i.e., classes
     * implementing this interface will in general want to derive from both
     * this Interface class as well as from the SimulatorAccess class.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class Interface : public Plugins::InterfaceBase
    {
      public:
        /**
         * Execute this postprocessor. Derived classes will implement this
         * function to do whatever they want to do to evaluate the solution at
         * the current time step.
         *
         * @param[in,out] statistics An object that contains statistics that
         * are collected throughout the simulation and that will be written to
         * an output file at the end of each time step. Postprocessors may
         * deposit data in these tables for later visualization or further
         * processing.
         *
         * @return A pair of strings that will be printed to the screen after
         * running the postprocessor in two columns; typically the first
         * column contains a description of what the data is and the second
         * contains a numerical value of this data. If there is nothing to
         * print, simply return two empty strings.
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) = 0;

        /**
         * A function that is used to indicate to the postprocessor manager which
         * other postprocessor(s) the current one depends upon. The returned
         * list contains the names (as strings, as you would write them in
         * the input file) of the postprocessors it requires. The manager
         * will ensure that these postprocessors are indeed used, even if
         * they were not explicitly listed in the input file, and are indeed
         * run <i>before</i> this postprocessor every time they are executed.
         *
         * The default implementation of this function returns an empty list.
         *
         * @note This mechanism is intended to support postprocessors that
         * do not want to recompute information that other postprocessors
         * already compute (or could compute, if they were run -- which
         * we can ensure using the current function). To do so, a postprocessor
         * of course needs to be able to access these other postprocessors.
         * This can be done by deriving your postprocessor from
         * SimulatorAccess, and then using the
         * SimulatorAccess::get_postprocess_manager() function, followed
         * by asking the resulting object via get_matching_active_plugin()
         * for a specific postprocessor object.
         */
        virtual
        std::list<std::string>
        required_other_postprocessors () const;
    };



    /**
     * A class that manages all objects that provide functionality to
     * postprocess solutions. It declares run time parameters for input files,
     * reads their values from such an input file, manages a list of all
     * postprocessors selected in the input file, and upon request through the
     * execute() function calls them in turn.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class Manager : public Plugins::ManagerBase<Interface<dim>>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Execute all of the postprocessor objects that have been requested
         * in the input file. These objects also fill the contents of the
         * statistics object.
         *
         * The function returns a concatenation of the text returned by the
         * individual postprocessors.
         */
        std::list<std::pair<std::string,std::string>>
        execute (TableHandler &statistics);

        /**
         * Go through the list of all postprocessors that have been selected
         * in the input file (and are consequently currently active) and return
         * true if one of them has the desired type specified by the template
         * argument.
         *
         * This function can only be called if the given template type (the first template
         * argument) is a class derived from the Interface class in this namespace.
         *
         * @deprecated Instead of this function, use the
         *   Plugins::ManagerBase::has_matching_active_plugin() and
         *   Plugins::ManagerBase::get_matching_active_plugin() functions of the base
         *   class of the current class.
         */
        template <typename PostprocessorType,
                  typename = typename std::enable_if_t<std::is_base_of<Interface<dim>,PostprocessorType>::value>>
        DEAL_II_DEPRECATED
        bool
        has_matching_postprocessor () const;

        /**
         * Go through the list of all postprocessors that have been selected
         * in the input file (and are consequently currently active) and see
         * if one of them has the type specified by the template
         * argument or can be cast to that type. If so, return a reference
         * to it. If no postprocessor is active that matches the given type,
         * throw an exception.
         *
         * This function can only be called if the given template type (the first template
         * argument) is a class derived from the Interface class in this namespace.
         *
         * @deprecated Instead of this function, use the
         *   Plugins::ManagerBase::has_matching_active_plugin() and
         *   Plugins::ManagerBase::get_matching_active_plugin() functions of the base
         *   class of the current class.
         */
        template <typename PostprocessorType,
                  typename = typename std::enable_if_t<std::is_base_of<Interface<dim>,PostprocessorType>::value>>
        DEAL_II_DEPRECATED
        const PostprocessorType &
        get_matching_postprocessor () const;

        /**
         * Declare the parameters of all known postprocessors, as well as of
         * ones this class has itself.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * This determines which postprocessor objects will be created; then
         * let these objects read their parameters as well.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * A function that is used to register postprocessor objects in such a
         * way that the Manager can deal with all of them without having to
         * know them by name. This allows the files in which individual
         * postprocessors are implement to register these postprocessors,
         * rather than also having to modify the Manager class by adding the
         * new postprocessor class.
         *
         * @param name The name under which this postprocessor is to be called
         * in parameter files.
         * @param description A text description of what this model does and
         * that will be listed in the documentation of the parameter file.
         * @param declare_parameters_function A pointer to a function that
         * declares the parameters for this postprocessor.
         * @param factory_function A pointer to a function that creates such a
         * postprocessor object and returns a pointer to it.
         */
        static
        void
        register_postprocessor (const std::string &name,
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
        DeclException1 (ExcPostprocessorNameNotFound,
                        std::string,
                        << "Could not find entry <"
                        << arg1
                        << "> among the names of registered postprocessors.");
    };


    /* -------------------------- inline and template functions ---------------------- */

    template <int dim>
    template <typename PostprocessorType, typename>
    inline
    bool
    Manager<dim>::has_matching_postprocessor () const
    {
      return this->template has_matching_active_plugin<PostprocessorType>();
    }



    template <int dim>
    template <typename PostprocessorType, typename>
    inline
    const PostprocessorType &
    Manager<dim>::get_matching_postprocessor () const
    {
      return this->template get_matching_active_plugin<PostprocessorType>();
    }


    /**
     * Given a class name, a name, and a description for the parameter file
     * for a postprocessor, register it with the aspect::Postprocess::Manager
     * class.
     *
     * @ingroup Postprocessing
     */
#define ASPECT_REGISTER_POSTPROCESSOR(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_POSTPROCESSOR_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::Postprocess::Interface<2>,classname<2>> \
    dummy_ ## classname ## _2d (&aspect::Postprocess::Manager<2>::register_postprocessor, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::Postprocess::Interface<3>,classname<3>> \
    dummy_ ## classname ## _3d (&aspect::Postprocess::Manager<3>::register_postprocessor, \
                                name, description); \
  }
  }
}


#endif
