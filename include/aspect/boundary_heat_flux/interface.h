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


#ifndef _aspect_boundary_heat_flux_interface_h
#define _aspect_boundary_heat_flux_interface_h

#include <aspect/plugins.h>
#include <aspect/material_model/interface.h>

namespace aspect
{
  /**
   * A namespace for the definition of things that have to do with describing
   * the boundary heat flux values.
   *
   * @ingroup BoundaryHeatFlux
   */
  namespace BoundaryHeatFlux
  {
    /**
     * Base class
     *
     * @ingroup BoundaryHeatFlux
     */
    template <int dim>
    class Interface : public Plugins::InterfaceBase
    {
      public:
        /**
         * Compute the heat flux for a list of quadrature points.
         *
         * The return value would typically be computed as the product of the thermal
         * conductivity @p material_model_outputs.thermal_conductivities[q] and the
         * temperature gradient at the boundary.
         *
         * @param boundary_indicator The boundary indicator of the part of the
         * boundary of the domain on which the point is located at which we
         * are requesting the fluid pressure gradients.
         * @param material_model_inputs The material property inputs.
         * @param material_model_outputs The material property outputs.
         * @param normal_vectors The normal vector at each quadrature point.
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
        std::vector<Tensor<1,dim>>
        heat_flux (const types::boundary_id boundary_indicator,
                   const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                   const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                   const std::vector<Tensor<1,dim>> &normal_vectors) const = 0;
    };


    /**
     * Register a boundary heat flux model so that it can be selected from
     * the parameter file.
     *
     * @param name A string that identifies the fluid pressure boundary model
     * @param description A text description of what this model does and that
     * will be listed in the documentation of the parameter file.
     * @param declare_parameters_function A pointer to a function that can be
     * used to declare the parameters that this boundary heat flux model
     * wants to read from input files.
     * @param factory_function A pointer to a function that can create an
     * object of this boundary heat flux model.
     *
     * @ingroup BoundaryHeatFlux
     */
    template <int dim>
    void
    register_boundary_heat_flux (const std::string &name,
                                 const std::string &description,
                                 void (*declare_parameters_function) (ParameterHandler &),
                                 std::unique_ptr<Interface<dim>> (*factory_function) ());

    /**
     * A function that given the name of a model returns a pointer to an
     * object that describes it. Ownership of the pointer is transferred to
     * the caller.
     *
     * The model object returned is not yet initialized and has not
     * read its runtime parameters yet.
     *
     * @ingroup BoundaryHeatFlux
     */
    template <int dim>
    std::unique_ptr<Interface<dim>>
    create_boundary_heat_flux (ParameterHandler &prm);


    /**
     * Declare the runtime parameters of the registered boundary heat flux
     * models.
     *
     * @ingroup BoundaryHeatFlux
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
     * for a boundary heat flux model, register it with the functions that
     * can declare their parameters and create these objects.
     *
     * @ingroup BoundaryHeatFlux
     */
#define ASPECT_REGISTER_BOUNDARY_HEAT_FLUX_MODEL(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_BOUNDARY_HEAT_FLUX_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::BoundaryHeatFlux::Interface<2>,classname<2>> \
    dummy_ ## classname ## _2d (&aspect::BoundaryHeatFlux::register_boundary_heat_flux<2>, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::BoundaryHeatFlux::Interface<3>,classname<3>> \
    dummy_ ## classname ## _3d (&aspect::BoundaryHeatFlux::register_boundary_heat_flux<3>, \
                                name, description); \
  }
  }
}


#endif
