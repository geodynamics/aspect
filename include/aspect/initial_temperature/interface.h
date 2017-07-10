/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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


#ifndef _aspect_initial_temperature_interface_h
#define _aspect_initial_temperature_interface_h

#include <aspect/plugins.h>
#include <aspect/utilities.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  template <int dim> class SimulatorAccess;

  /**
   * A namespace in which we define everything that has to do with defining
   * the initial conditions.
   *
   * @ingroup InitialTemperatures
   */
  namespace InitialTemperature
  {
    using namespace dealii;

    /**
     * A base class for parameterizations of initial conditions.
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class Interface
    {
      public:
        /**
         * Destructor. Made virtual to enforce that derived classes also have
         * virtual destructors.
         */
        virtual ~Interface();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after
         * the SimulatorAccess (if applicable) is initialized.
         */
        virtual
        void
        initialize ();

        /**
         * Return the initial temperature as a function of position.
         */
        virtual
        double initial_temperature (const Point<dim> &position) const = 0;


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
     * A class that manages all initial condition objects.
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class Manager : public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Destructor. Made virtual since this class has virtual member
         * functions.
         */
        virtual ~Manager ();

        /**
         * Declare the parameters of all known initial conditions plugins, as
         * well as of ones this class has itself.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * This determines which initial conditions objects will be created; then
         * let these objects read their parameters as well.
         */
        void
        parse_parameters (ParameterHandler &prm);

        /**
         * A function that calls the initial_temperature functions of all the
         * individual initial condition objects and adds up the values of the
         * individual calls.
         */
        double
        initial_temperature (const Point<dim> &position) const;

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
        static
        void
        register_initial_temperature (const std::string &name,
                                      const std::string &description,
                                      void (*declare_parameters_function) (ParameterHandler &),
                                      Interface<dim> *(*factory_function) ());


        /**
         * Return a list of names of all initial temperature models currently
         * used in the computation, as specified in the input file.
         */
        const std::vector<std::string> &
        get_active_initial_temperature_names () const;

        /**
         * Return a list of pointers to all initial temperature models
         * currently used in the computation, as specified in the input file.
         */
        const std::list<std_cxx11::shared_ptr<Interface<dim> > > &
        get_active_initial_temperature_conditions () const;

        /**
         * Go through the list of all initial temperature models that have been selected in
         * the input file (and are consequently currently active) and see if one
         * of them has the desired type specified by the template argument. If so,
         * return a pointer to it. If no initial temperature model is active
         * that matches the given type, return a NULL pointer.
         */
        template <typename InitialTemperatureType>
        InitialTemperatureType *
        find_initial_temperature_model () const;

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
        DeclException1 (ExcInitialTemperatureNameNotFound,
                        std::string,
                        << "Could not find entry <"
                        << arg1
                        << "> among the names of registered initial temperature objects.");
      private:
        /**
         * A list of initial temperature objects that have been requested in the
         * parameter file.
         */
        std::list<std_cxx11::shared_ptr<Interface<dim> > > initial_temperature_objects;

        /**
         * A list of names of initial temperature objects that have been requested
         * in the parameter file.
         */
        std::vector<std::string> model_names;

        /**
         * A list of enums of initial temperature operators that have been
         * requested in the parameter file. Each entry is used to modify the
         * initial temperature field with the values from the associated plugin
         * in model_names.
         */
        std::vector<aspect::Utilities::Operator> model_operators;
    };



    template <int dim>
    template <typename InitialTemperatureType>
    inline
    InitialTemperatureType *
    Manager<dim>::find_initial_temperature_model () const
    {
      for (typename std::list<std_cxx11::shared_ptr<Interface<dim> > >::const_iterator
           p = initial_temperature_objects.begin();
           p != initial_temperature_objects.end(); ++p)
        if (InitialTemperatureType *x = dynamic_cast<InitialTemperatureType *> ( (*p).get()) )
          return x;
      return NULL;
    }



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
     * Given a class name, a name, and a description for the parameter file
     * for a initial conditions model, register it with the functions that can
     * declare their parameters and create these objects.
     *
     * @ingroup InitialTemperatures
     */
#define ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::InitialTemperature::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::InitialTemperature::Manager<2>::register_initial_temperature, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::InitialTemperature::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::InitialTemperature::Manager<3>::register_initial_temperature, \
                                name, description); \
  }
  }
}


#endif
