/*
  Copyright (C) 2013 by the authors of the ASPECT code.

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


#ifndef __aspect__adiabatic_conditions_interface_h
#define __aspect__adiabatic_conditions_interface_h

#include <aspect/plugins.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/distributed/tria.h>


namespace aspect
{
  /**
   * A namespace for the definition of things that have to do with describing
   * the calculation of a reference adiabatic profile.
   *
   * @ingroup AdiabaticConditions
   */
  namespace AdiabaticConditions
  {
    using namespace dealii;

    /**
     * Base class for classes that describe adiabatic conditions, i.e. that
     * starts at the top of the domain and integrate pressure and temperature
     * as we go down into depth. There are several ways to do this (time-
     * dependent or constant, using laterally averaged values or a reference
     * profile), therefore we allow for user written plugins.
     *
     * @ingroup AdiabaticConditions
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
         * Some plugins need to know whether the adiabatic conditions are
         * already calculated. Namely all plugins that are needed to create
         * the adiabatic conditions but themselves depedend on the adiabatic
         * profile. Utilizing this function they may behave differently on
         * initialization of the adiabatic conditions and at model runtime.
         */
        virtual
        bool
        is_initialized () const = 0;

        /**
         * Compute the adiabatic conditions along a vertical transect of the
         * geometry based on the given material model and other quantities.
         * This function is called at every new timestep.
         */
        virtual
        void update ();

        /**
         * Return the adiabatic temperature at a given point of the domain.
         */
        virtual
        double temperature (const Point<dim> &p) const = 0;

        /**
         * Return the adiabatic temperature profile as a vector of values
         * corresponding to increasing depth.
         *
         * @param values The output vector of depth averaged values. The
         * function takes the pre-existing size of this vector as the number
         * of depth slices.
         */
        virtual
        void get_adiabatic_temperature_profile(std::vector<double> &values) const = 0;

        /**
         * Return the adiabatic pressure at a given point of the domain.
         */
        virtual
        double pressure (const Point<dim> &p) const = 0;


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
     * Register a adiabatic conditions model so that it can be selected from
     * the parameter file.
     *
     * @param name A string that identifies the adiabatic conditions model
     * @param description A text description of what this model does and that
     * will be listed in the documentation of the parameter file.
     * @param declare_parameters_function A pointer to a function that can be
     * used to declare the parameters that this model wants to read from input
     * files.
     * @param factory_function A pointer to a function that can create an
     * object of this adiabatic conditions model.
     *
     * @ingroup AdiabaticConditions
     */
    template <int dim>
    void
    register_adiabatic_conditions (const std::string &name,
                                   const std::string &description,
                                   void (*declare_parameters_function) (ParameterHandler &),
                                   Interface<dim> *(*factory_function) ());

    /**
     * A function that given the name of a model returns a pointer to an
     * object that describes it. Ownership of the pointer is transferred to
     * the caller.
     *
     * The model object returned is not yet initialized and has not read its
     * runtime parameters yet.
     *
     * @ingroup AdiabaticConditions
     */
    template <int dim>
    Interface<dim> *
    create_adiabatic_conditions (ParameterHandler &prm);


    /**
     * Declare the runtime parameters of the registered adiabatic conditions
     * model.
     *
     * @ingroup AdiabaticConditions
     */
    template <int dim>
    void
    declare_parameters (ParameterHandler &prm);


    /**
     * Given a class name, a name, and a description for the parameter file
     * for a adiabatic conditions model, register it with the functions that
     * can declare their parameters and create these objects.
     *
     * @ingroup AdiabaticConditions
     */
#define ASPECT_REGISTER_ADIABATIC_CONDITIONS_MODEL(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_ADIABATIC_CONDITIONS_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::AdiabaticConditions::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::AdiabaticConditions::register_adiabatic_conditions<2>, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::AdiabaticConditions::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::AdiabaticConditions::register_adiabatic_conditions<3>, \
                                name, description); \
  }
  }
}


#endif
