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


#ifndef _aspect_boundary_fluid_pressure_interface_h
#define _aspect_boundary_fluid_pressure_interface_h

#include <aspect/plugins.h>
#include <aspect/material_model/interface.h>

namespace aspect
{
  /**
   * A namespace for the definition of things that have to do with describing
   * the boundary values for fluid pressure for computations with melt
   * transport.
   *
   * @ingroup BoundaryFluidPressures
   */
  namespace BoundaryFluidPressure
  {
    /**
     * Base class
     *
     * @ingroup BoundaryFluidPressures
     */
    template <int dim>
    class Interface : public Plugins::InterfaceBase
    {
      public:
        /**
         * Compute the component of the gradient of the fluid pressure
         * in the direction normal to a boundary for a list of quadrature
         * points.
         *
         * The return value can typically contain @p material_model_outputs.fluid_densities[q]
         * or @p material_model_outputs.densities[q], multiplied by the gravity vector
         * and dotted with the normal.
         * If the solid density is used, fluid is only flowing in or out due to differences in
         * dynamic pressure, if the fluid density is used, melt flows in with the same velocity
         * as inflowing solid material.
         *
         * @param boundary_indicator The boundary indicator of the part of the
         * boundary of the domain on which the point is located at which we
         * are requesting the fluid pressure gradients.
         * @param material_model_inputs The material property inputs.
         * @param material_model_outputs The material property outputs.
         * @param normal_vectors A normal vector for each point.
         * @param fluid_pressure_gradient_outputs Result to be filled.
         */
        virtual
        void fluid_pressure_gradient (
          const types::boundary_id boundary_indicator,
          const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
          const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
          const std::vector<Tensor<1,dim>> &normal_vectors,
          std::vector<double> &fluid_pressure_gradient_outputs
        ) const = 0;
    };


    /**
     * Register a fluid pressure boundary model so that it can be selected from
     * the parameter file.
     *
     * @param name A string that identifies the fluid pressure boundary model
     * @param description A text description of what this model does and that
     * will be listed in the documentation of the parameter file.
     * @param declare_parameters_function A pointer to a function that can be
     * used to declare the parameters that this fluid pressure boundary model
     * wants to read from input files.
     * @param factory_function A pointer to a function that can create an
     * object of this fluid pressure boundary model.
     *
     * @ingroup BoundaryFluidPressures
     */
    template <int dim>
    void
    register_boundary_fluid_pressure (const std::string &name,
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
     * @ingroup BoundaryFluidPressures
     */
    template <int dim>
    std::unique_ptr<Interface<dim>>
    create_boundary_fluid_pressure (ParameterHandler &prm);


    /**
     * Declare the runtime parameters of the registered fluid pressure boundary
     * models.
     *
     * @ingroup BoundaryFluidPressures
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
     * for a fluid pressure boundary model, register it with the functions that
     * can declare their parameters and create these objects.
     *
     * @ingroup BoundaryFluidPressures
     */
#define ASPECT_REGISTER_BOUNDARY_FLUID_PRESSURE_MODEL(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_BOUNDARY_FLUID_PRESSURE_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::BoundaryFluidPressure::Interface<2>,classname<2>> \
    dummy_ ## classname ## _2d (&aspect::BoundaryFluidPressure::register_boundary_fluid_pressure<2>, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::BoundaryFluidPressure::Interface<3>,classname<3>> \
    dummy_ ## classname ## _3d (&aspect::BoundaryFluidPressure::register_boundary_fluid_pressure<3>, \
                                name, description); \
  }
  }
}


#endif
