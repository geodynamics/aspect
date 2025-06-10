/*
  Copyright (C) 2015 - 2024 by the authors of the ASPECT code.

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


#ifndef _aspect_boundary_convective_heating_interface_h
#define _aspect_boundary_convective_heating_interface_h

#include <aspect/plugins.h>
#include <aspect/material_model/interface.h>
#include <aspect/boundary_temperature/interface.h>
#include <aspect/boundary_heat_flux/interface.h>
#include <aspect/utilities.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  template <int dim> class SimulatorAccess;

  /**
   * A namespace for the definition of things that have to do with describing
   * the boundary heat flux values.
   *
   * @ingroup BoundaryConvectiveHeating
   */
  namespace BoundaryConvectiveHeating
  {
    /**
     * Base class
     *
     * @ingroup BoundaryConvectiveHeating
     */
    template <int dim>
    class Interface : public Plugins::InterfaceBase
    {
      public:
        /**
         * Compute the heat transfer coefficient for a list of quadrature points.
         *
         * TODO
         *
         * @param boundary_indicator The boundary indicator of the part of the
         * boundary of the domain on which the point is located at which we
         * are requesting the fluid pressure gradients.
         * @param material_model_inputs The material property inputs.
         * @param material_model_outputs The material property outputs.
         *
         * @return A vector of heat flux vectors at the evaluation points.
         *   For historical reasons, the function is asked
         *   to provide the heat flux as a vector, even though the place where the
         *   heat flux is used only uses the component of this vector that is
         *   to the boundary (which it computes by taking the dot product *normal*
         *   between the returned vector and the normal vector). Because there are
         *   situations where all you can do is compute the normal heat flux as a
         *   scalar, the `heat_flux()` function also receives the normal vector as
         *   an input argument. As a consequence, one way for the function to
         *   compute the required heat flux vector is to compute the scalar heat
         *   flux and multiply it by the normal vector.
         */
        virtual
        std::vector<double>
        heat_transfer_coefficient (const types::boundary_id boundary_indicator,
                                   const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                                   const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs) const = 0;
    };



    /**
     * A class that manages all boundary convective heating objects.
     *
     * @ingroup BoundaryConvectiveHeating
     */
    template <int dim>
    class Manager : public Plugins::ManagerBase<Interface<dim>>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Declare the parameters of all known boundary convective heating plugins, as
         * well as the ones this class has itself.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * This determines which boundary convective heating, boundary temperature
         * and boundary heat flux objects will be created;
         * then let these objects read their parameters as well.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * A function that calls the boundary_temperature functions of all the
         * individual boundary temperature objects and uses the stored operators
         * to combine them.
         */
        double
        boundary_temperature (const types::boundary_id boundary_indicator,
                              const Point<dim> &position) const;

        /**
         * A function that calls the heat_flux functions of the
         * individual boundary heat flux object.
         */
        std::vector<Tensor<1,dim>>
        heat_flux (const types::boundary_id boundary_indicator,
                   const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                   const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                   const std::vector<Tensor<1,dim>> &normal_vectors) const;

        /**
         * A function that calls the heat_transfer_coefficient functions of the
         * individual boundary convection heating object.
         */
        std::vector<double>
        heat_transfer_coefficient (const types::boundary_id boundary_indicator,
                                   const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                                   const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs) const;

        /**
         * A function that is used to register boundary convective heating objects in such
         * a way that the Manager can deal with all of them without having to
         * know them by name. This allows the files in which individual
         * plugins are implemented to register these plugins, rather than also
         * having to modify the Manager class by adding the new boundary
         * convective heating plugin class.
         *
         * @param name A string that identifies the boundary temperature model
         * @param description A text description of what this model does and that
         * will be listed in the documentation of the parameter file.
         * @param declare_parameters_function A pointer to a function that can be
         * used to declare the parameters that this boundary temperature model
         * wants to read from input files.
         * @param factory_function A pointer to a function that can create an
         * object of this boundary temperature model.
         */
        static
        void
        register_boundary_convective_heating (const std::string &name,
                                              const std::string &description,
                                              void (*declare_parameters_function) (ParameterHandler &),
                                              std::unique_ptr<Interface<dim>> (*factory_function) ());

        static
        void
        register_boundary_temperature (const std::string &name,
                                       const std::string &description,
                                       void (*declare_parameters_function) (ParameterHandler &),
                                       std::unique_ptr<BoundaryTemperature::Interface<dim>> (*factory_function) ());

        static
        void
        register_boundary_heat_flux (const std::string &name,
                                     const std::string &description,
                                     void (*declare_parameters_function) (ParameterHandler &),
                                     std::unique_ptr<BoundaryHeatFlux::Interface<dim>> (*factory_function) ());

        /*
         * Return a set of boundary indicators for which convective heating
         * is prescribed.
         */
        const std::set<types::boundary_id> &
        get_convective_heating_boundary_indicators() const;


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
        DeclException1 (ExcBoundaryTemperatureNameNotFound,
                        std::string,
                        << "Could not find entry <"
                        << arg1
                        << "> among the names of registered boundary temperature objects.");
      private:
        /**
         * A list of enums of boundary temperature operators that have been
         * requested in the parameter file. Each name is associated
         * with a model_name, and is used to modify the temperature
         * boundary with the values from the current plugin.
         */
        std::vector<aspect::Utilities::Operator> model_operators;

        /**
         * A set of boundary ids on which the boundary_temperature_objects
         * will be applied.
         */
        std::set<types::boundary_id> convective_heating_boundary_indicators;

        /**
        * Names of temperature plugins.
        */
        std::vector<std::string> temperature_plugin_names;

        /**
         * Names of heat flux plugins.
         */
        std::vector<std::string> heat_flux_plugin_names;

        /**
         * Objects for temperature plugins.
         */
        std::list<std::unique_ptr<aspect::BoundaryTemperature::Interface<dim>>> temperature_plugin_objects;

        /**
         * Objects for heat flux plugins.
         */
        std::list<std::unique_ptr<aspect::BoundaryHeatFlux::Interface<dim>>> heat_flux_plugin_objects;
    };


    /**
    * Given a class name, a name, and a description for the parameter file
    * for a boundary convective heating model, register it with the functions that
    * can declare their parameters and create these objects.
    *
    * @ingroup BoundaryConvectiveHeating
    */
#define ASPECT_REGISTER_BOUNDARY_CONVECTIVE_HEATING_MODEL(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_BOUNDARY_CONVECTIVE_HEATING_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::BoundaryConvectiveHeating::Interface<2>,classname<2>> \
    dummy_ ## classname ## _2d (&aspect::BoundaryConvectiveHeating::Manager<2>::register_boundary_convective_heating, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::BoundaryConvectiveHeating::Interface<3>,classname<3>> \
    dummy_ ## classname ## _3d (&aspect::BoundaryConvectiveHeating::Manager<3>::register_boundary_convective_heating, \
                                name, description); \
  }
  }
}


#endif
