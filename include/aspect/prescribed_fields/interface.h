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


#ifndef _aspect_prescribed_fields_interface_h
#define _aspect_prescribed_fields_interface_h

#include <aspect/plugins.h>
#include <aspect/utilities.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/function.h>

#include <boost/core/demangle.hpp>
#include <typeinfo>


namespace aspect
{
  template <int dim> class SimulatorAccess;

  /**
   * A namespace in which we define everything that has to do with defining
   * the initial conditions.
   *
   * @ingroup PrescribedFieldss
   */
  namespace PrescribedFields
  {
    /**
    * This plugin allows the user to prescribe internal fields and can be
    * thought of as temperature, velocities, etc. equivalent of the initial
    * conditions plugin.
    *
    *
    * @ingroup PrescribedFields
    */
    template <int dim>
    class Interface : public Plugins::InterfaceBase
    {
      public:

        /**
         * A function that is called at the beginning of each time step to
         * indicate what the model time is for which the prescribed fields
         * are going to be evaluated.
         */
        virtual
        void update(const SimulatorAccess<dim> &simulator_access) = 0;

        /**
         * By having access to the simulator and current constraints, add new constraints
         * based on individual plugins
         */
        virtual
        void constrain_internal_fields (const SimulatorAccess<dim> &simulator_access,
                                        AffineConstraints<double> &current_constraints) = 0;
    };

    /**
    * @ingroup PrescribedFields
    */
    /**
    * A function that is used to register initial temperature objects in such
    * a way that the Manager can deal with all of them without having to
    * know them by name. This allows the files in which individual
    * plugins are implemented to register these plugins, rather than also
    * having to modify the Manager class by adding the new initial
    * temperature plugin class.
    *
    * @param name A string that identifies the initial temperature model
    * @param description A text description of what this model does and that
    * will be listed in the documentation of the parameter file.
    * @param declare_parameters_function A pointer to a function that can be
    * used to declare the parameters that this initial temperature model
    * wants to read from input files.
    * @param factory_function A pointer to a function that can create an
    * object of this initial temperature model.
    */
    template <int dim>
    void
    register_prescribed_fields_model (const std::string &name,
                                      const std::string &description,
                                      void (*declare_parameters_function) (ParameterHandler &),
                                      std::unique_ptr<Interface<dim>> (*factory_function) ());


    /**
     * Declare the parameters of all known initial conditions plugins, as
     * well as of ones this class has itself.
     */
    template <int dim>
    void
    declare_parameters(ParameterHandler &prm);

    /**
     * A wrapper of the declaration that handles the dimension explicitly.
     * This format is needed for connecting with SimulatorSignals
     */
    void
    declare_parameters_signal(const unsigned int dim,
                              ParameterHandler &prm);


    /**
     * Read the parameters this class declares from the parameter file.
     * This determines which initial conditions objects will be created; then
     * let these objects read their parameters as well.
     */
    template <int dim>
    void parse_parameters(const Parameters<dim> &,
                          ParameterHandler &prm);


    /**
     * For the current plugin subsystem, write a connection graph of all of the
     * plugins we know about, in the format that the
     * programs dot and neato understand. This allows for a visualization of
     * how all of the plugins that ASPECT knows about are interconnected, and
     * connect to other parts of the ASPECT code.
     *
     * @param output_stream The stream to write the output to.
     */
    void
    write_plugin_graph (std::ostream &output_stream);

    /**
     * Exception.
     */
    DeclException1 (ExcPrescribedFieldsNameNotFound,
                    std::string,
                    << "Could not find entry <"
                    << arg1
                    << "> among the names of registered initial temperature objects.");

    /**
     * whether internal fiels are prescribed
    */
    bool prescribe_internal_fields;

    /**
     * A list of names used in the input file to identify plugins,
     * corresponding to the plugin objects stored in the previous variable.
     */
    std::vector<std::string> plugin_names;

    /**
     * Objects for prescribed fields plugins.
     * Because we don't know what dimension we are handing,
     * we keep both 2d and 3d instances.
     */
    std::list<std::unique_ptr<aspect::PrescribedFields::Interface<2>>> plugin_objects_2d;
    std::list<std::unique_ptr<aspect::PrescribedFields::Interface<3>>> plugin_objects_3d;


    /**
     * Return a string that consists of the names of initial temperature models that can
     * be selected. These names are separated by a vertical line '|' so
     * that the string can be an input to the deal.II classes
     * Patterns::Selection or Patterns::MultipleSelection.
     */
    template <int dim>
    std::string
    get_valid_model_names_pattern ();

    /**
     * A function that calls the constrain_internal_fields functions of all the
     * individual prescribed fields objects and assign new constrains.
     */
    template <int dim>
    void constrain_all_prescribed_internal_fields (const SimulatorAccess<dim> &simulator_access,
                                                   AffineConstraints<double> &current_constraints);
    /**
    * Given a class name, a name, and a description for the parameter file
    * for a prescribed Stokes solution model, register it with the functions
    * that can declare their parameters and create these objects.
    *
    * @ingroup PrescribedFields
    */
#define ASPECT_REGISTER_PRESCRIBED_FIELDS(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_PRESCRIBED_FIELDS_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::PrescribedFields::Interface<2>,classname<2>> \
    dummy_ ## classname ## _2d (&aspect::PrescribedFields::register_prescribed_fields_model<2>, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::PrescribedFields::Interface<3>,classname<3>> \
    dummy_ ## classname ## _3d (&aspect::PrescribedFields::register_prescribed_fields_model<3>, \
                                name, description); \
  }
  }
}


#endif
