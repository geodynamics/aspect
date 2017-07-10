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


#ifndef _aspect_initial_composition_interface_h
#define _aspect_initial_composition_interface_h

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
   * the initial composition.
   *
   * @ingroup InitialCompositions
   */
  namespace InitialComposition
  {
    using namespace dealii;

    /**
     * A base class for parameterizations of initial conditions.
     *
     * @ingroup InitialCompositions
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
         * Return the initial composition as a function of position.
         */
        virtual
        double initial_composition (const Point<dim> &position, const unsigned int n_comp) const = 0;


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
     * A class that manages all initial composition objects.
     *
     * @ingroup InitialCompositions
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
         * Declare the parameters of all known initial composition plugins, as
         * well as of ones this class has itself.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * This determines which initial composition objects will be created;
         * then let these objects read their parameters as well.
         */
        void
        parse_parameters (ParameterHandler &prm);

        /**
         * A function that calls the initial_composition functions of all
         * individual initial composition objects for the composition given
         * by @p n_comp, and adds up the values of the
         * individual calls.
         */
        double
        initial_composition (const Point<dim> &position,
                             const unsigned int n_comp) const;

        /**
         * A function that is used to register initial composition objects in
         * such a way that the Manager can deal with all of them without having
         * to know them by name. This allows the files in which individual
         * plugins are implemented to register these plugins, rather than also
         * having to modify the Manager class by adding the new initial
         * composition plugin class.
         *
         * @param name A string that identifies the initial composition model
         * @param description A text description of what this model does and that
         * will be listed in the documentation of the parameter file.
         * @param declare_parameters_function A pointer to a function that can be
         * used to declare the parameters that this initial composition model
         * wants to read from input files.
         * @param factory_function A pointer to a function that can create an
         * object of this initial composition model.
         */
        static
        void
        register_initial_composition (const std::string &name,
                                      const std::string &description,
                                      void (*declare_parameters_function) (ParameterHandler &),
                                      Interface<dim> *(*factory_function) ());


        /**
         * Return a list of names of all initial composition models currently
         * used in the computation, as specified in the input file.
         */
        const std::vector<std::string> &
        get_active_initial_composition_names () const;

        /**
         * Return a list of pointers to all initial composition models
         * currently used in the computation, as specified in the input file.
         */
        const std::list<std_cxx11::shared_ptr<Interface<dim> > > &
        get_active_initial_composition_conditions () const;

        /**
         * Go through the list of all initial composition models that have been selected in
         * the input file (and are consequently currently active) and see if one
         * of them has the desired type specified by the template argument. If so,
         * return a pointer to it. If no initial composition model is active that matches the
         * given type, return a NULL pointer.
         */
        template <typename InitialCompositionType>
        InitialCompositionType *
        find_initial_composition_model () const;

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
        DeclException1 (ExcInitialCompositionNameNotFound,
                        std::string,
                        << "Could not find entry <"
                        << arg1
                        << "> among the names of registered initial composition objects.");
      private:
        /**
         * A list of initial composition objects that have been requested in the
         * parameter file.
         */
        std::list<std_cxx11::shared_ptr<Interface<dim> > > initial_composition_objects;

        /**
         * A list of names of initial composition objects that have been requested
         * in the parameter file.
         */
        std::vector<std::string> model_names;

        /**
         * A list of enums of initial composition operators that have been
         * requested in the parameter file. Each entry is used to modify the
         * initial compositional field with the values from the associated plugin
         * in model_names.
         */
        std::vector<aspect::Utilities::Operator> model_operators;
    };



    /**
     * A function that given the name of a model returns a pointer to an
     * object that describes it. Ownership of the pointer is transferred to
     * the caller.
     *
     * The model object returned is not yet initialized and has not read its
     * runtime parameters yet.
     *
     * @ingroup InitialCompositions
     */
    template <int dim>
    template <typename InitialCompositionType>
    inline
    InitialCompositionType *
    Manager<dim>::find_initial_composition_model () const
    {
      for (typename std::list<std_cxx11::shared_ptr<Interface<dim> > >::const_iterator
           p = initial_composition_objects.begin();
           p != initial_composition_objects.end(); ++p)
        if (InitialCompositionType *x = dynamic_cast<InitialCompositionType *> ( (*p).get()) )
          return x;
      return NULL;
    }



    /**
     * Return a string that consists of the names of initial composition models that can
     * be selected. These names are separated by a vertical line '|' so
     * that the string can be an input to the deal.II classes
     * Patterns::Selection or Patterns::MultipleSelection.
     *
     * @ingroup InitialCompositions
     */
    template <int dim>
    std::string
    get_valid_model_names_pattern ();



    /**
     * Given a class name, a name, and a description for the parameter file
     * for a initial composition model, register it with the functions that can
     * declare their parameters and create these objects.
     *
     * @ingroup InitialCompositions
     */
#define ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::InitialComposition::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::InitialComposition::Manager<2>::register_initial_composition, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::InitialComposition::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::InitialComposition::Manager<3>::register_initial_composition, \
                                name, description); \
  }
  }
}


#endif
