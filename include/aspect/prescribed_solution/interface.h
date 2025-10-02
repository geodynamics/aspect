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


#ifndef _aspect_prescribed_solution_interface_h
#define _aspect_prescribed_solution_interface_h

#include <aspect/plugins.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>

#include <boost/core/demangle.hpp>
#include <typeinfo>


namespace aspect
{
  /**
   * A namespace in which we define everything that has to do with prescribing
   * the solution in certain regions.
   *
   * @ingroup PrescribedSolution
   */
  namespace PrescribedSolution
  {
    /**
    * This plugin allows the user to prescribe fields and can be
    * thought of as temperature, velocities, etc. equivalent of the initial
    * conditions plugin.
    *
    * @ingroup PrescribedSolution
    */
    template <int dim>
    class Interface : public Plugins::InterfaceBase
    {
      public:
        /**
         * A function that walks through all active DoF @p positions in a @p cell and determines
         * whether each of them @p should_be_constrained based on the @p component_indices.
         * If a constraint is required, its value is computed and stored in the @p solution.
         */
        virtual
        void constrain_solution (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                 const std::vector<Point<dim>> &positions,
                                 const std::vector<unsigned int> &component_indices,
                                 std::vector<bool> &should_be_constrained,
                                 std::vector<double> &solution) = 0;
    };

    /**
     * A class that manages all prescribed solution objects.
     *
     * @ingroup PrescribedSolution
     */
    template <int dim>
    class Manager : public Plugins::ManagerBase<Interface<dim>>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Exception.
         */
        DeclException1 (ExcPrescribedSolutionNameNotFound,
                        std::string,
                        << "Could not find entry <"
                        << arg1
                        << "> among the names of registered prescribed solution objects.");

        /**
         * A function that calls the constrain_solution functions of all the
         * individual plugin objects and integrates their return values into
         * the @p current_constraints object.
         */
        void constrain_solution (AffineConstraints<double> &current_constraints) const;

        /**
         * A function that is used to register prescribed solution objects in such
         * a way that the Manager can deal with all of them without having to
         * know them by name. This allows the files in which individual
         * plugins are implemented to register these plugins, rather than also
         * having to modify the Manager class by adding the new prescribed
         * solution plugin class.
         *
         * @param name A string that identifies the prescribed solution model
         * @param description A text description of what this model does and that
         * will be listed in the documentation of the parameter file.
         * @param declare_parameters_function A pointer to a function that can be
         * used to declare the parameters that this prescribed solution model
         * wants to read from input files.
         * @param factory_function A pointer to a function that can create an
         * object of this prescribed solution model.
         */
        static
        void
        register_prescribed_solution_plugin (const std::string &name,
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
        * Declare the parameters of all known prescribed solution plugins, as
        * well as the ones this class has itself.
        */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * This determines which prescribed solution objects will be created;
         * then let these objects read their parameters as well.
         */
        void
        parse_parameters (ParameterHandler &prm) override;
    };


    /**
     * Return a string that consists of the names of prescribed solution models that can
     * be selected. These names are separated by a vertical line '|' so
     * that the string can be an input to the deal.II classes
     * Patterns::Selection or Patterns::MultipleSelection.
     */
    template <int dim>
    std::string
    get_valid_model_names_pattern ();

    /**
    * Given a class name, a name, and a description for the parameter file
    * for a prescribed solution model, register it with the functions
    * that can declare their parameters and create these objects.
    *
    * @ingroup PrescribedSolution
    */
#define ASPECT_REGISTER_PRESCRIBED_SOLUTION(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_PRESCRIBED_SOLUTION_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::PrescribedSolution::Interface<2>,classname<2>> \
    dummy_ ## classname ## _2d (&aspect::PrescribedSolution::Manager<2>::register_prescribed_solution_plugin, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::PrescribedSolution::Interface<3>,classname<3>> \
    dummy_ ## classname ## _3d (&aspect::PrescribedSolution::Manager<3>::register_prescribed_solution_plugin, \
                                name, description); \
  }
  }
}


#endif
