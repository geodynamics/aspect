/*
  Copyright (C) 2013 - 2024 by the authors of the ASPECT code.

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


#ifndef _aspect_adiabatic_conditions_interface_h
#define _aspect_adiabatic_conditions_interface_h

#include <aspect/plugins.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/simulator_access.h>
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
    /**
     * Base class for classes that describe adiabatic conditions,
     * i.e. that starts at the top of the domain and integrate
     * pressure and temperature as we go down into depth. There are
     * several ways to do this (time-dependent or constant, using
     * laterally averaged values or a reference profile), therefore we
     * allow for user written plugins.
     *
     * @ingroup AdiabaticConditions
     */
    template <int dim>
    class Interface: public SimulatorAccess<dim>, public Plugins::InterfaceBase
    {
      public:
        /**
         * Some plugins need to know whether the adiabatic conditions are
         * already calculated. Namely all plugins that are needed to create
         * the adiabatic conditions but themselves depend on the adiabatic
         * profile. Utilizing this function they may behave differently on
         * initialization of the adiabatic conditions and at model runtime.
         */
        virtual
        bool
        is_initialized () const = 0;

        /**
         * Return the adiabatic temperature at a given point of the domain.
         */
        virtual
        double temperature (const Point<dim> &p) const = 0;

        /**
         * Return the adiabatic pressure at a given point of the domain.
         */
        virtual
        double pressure (const Point<dim> &p) const = 0;

        /**
         * Return the reference_density at a given point of the domain.
         */
        virtual
        double density (const Point<dim> &p) const = 0;

        /**
         * Return the derivative of the density with respect to depth
         * at the given point @p p.
         */
        virtual
        double density_derivative (const Point<dim> &p) const = 0;

        /**
         * Return the adiabatic temperature profile as a vector of values
         * corresponding to increasing depth.
         *
         * @param values The output vector of depth averaged values. The
         * function takes the pre-existing size of this vector as the number
         * of depth slices.
         *
         * @deprecated: This function is deprecated.
         * Use the function temperature() for specific positions instead.
         */
        DEAL_II_DEPRECATED
        virtual
        void get_adiabatic_temperature_profile(std::vector<double> &values) const;

        /**
         * Like get_adiabatic_temperature_profile() but for the pressure.
         *
         * @deprecated: This function is deprecated.
         * Use the function pressure() for specific positions instead.
         */
        DEAL_II_DEPRECATED
        virtual
        void get_adiabatic_pressure_profile(std::vector<double> &values) const;

        /**
         * Like get_adiabatic_temperature_profile() but for the density.
         *
         * @deprecated: This function is deprecated.
         * Use the function density() for specific positions instead.
         */
        DEAL_II_DEPRECATED
        virtual
        void get_adiabatic_density_profile(std::vector<double> &values) const;

        /**
         * Like get_adiabatic_temperature_profile() but for the density derivative.
         *
         * @deprecated: This function is deprecated.
         * Use the function density_derivative() for specific positions instead.
         */
        DEAL_II_DEPRECATED
        virtual
        void get_adiabatic_density_derivative_profile(std::vector<double> &values) const;
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
                                   std::unique_ptr<Interface<dim>> (*factory_function) ());

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
    std::unique_ptr<Interface<dim>>
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
     * For the current plugin subsystem, write a connection graph of all of the
     * plugins we know about, in the format that the
     * programs dot and neato understand. This allows for a visualization of
     * how all of the plugins that ASPECT knows about are interconnected, and
     * connect to other parts of the ASPECT code.
     *
     * @param output_stream The stream to write the output to.
     */
    template <int dim>
    void
    write_plugin_graph (std::ostream &output_stream);


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
    aspect::internal::Plugins::RegisterHelper<aspect::AdiabaticConditions::Interface<2>,classname<2>> \
    dummy_ ## classname ## _2d (&aspect::AdiabaticConditions::register_adiabatic_conditions<2>, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::AdiabaticConditions::Interface<3>,classname<3>> \
    dummy_ ## classname ## _3d (&aspect::AdiabaticConditions::register_adiabatic_conditions<3>, \
                                name, description); \
  }
  }
}


#endif
