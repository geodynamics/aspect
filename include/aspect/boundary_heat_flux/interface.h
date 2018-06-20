/*
  Copyright (C) 2015 - 2017 by the authors of the ASPECT code.

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
   * the boundary values for fluid pressure for computations with melt
   * transport.
   *
   * @ingroup BoundaryHeatFlux
   */
  namespace BoundaryHeatFlux
  {
    using namespace dealii;

    /**
     * Base class
     *
     * @ingroup BoundaryHeatFlux
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
         * beginning of the program after parse_parameters is run and after the
         * SimulatorAccess (if applicable) is initialized.
         */
        virtual void initialize ();

        /**
         * A function that is called at the beginning of each time step. The
         * default implementation of the function does nothing, but derived
         * classes that need more elaborate setups for a given time step may
         * overload the function.
         *
         * The point of this function is to allow complex boundary heat flux
         * models to do an initialization step once at the beginning of each
         * time step. An example would be a model that needs to call an
         * external program to compute heat flux change at bottom.
         */
        virtual
        void
        update ();

        /**
         * Compute the component of the gradient of the fluid pressure
         * in the direction normal to a boundary for a list of quadrature
         * points.
         *
         * The return value can typically contain @p material_model_outputs.fluid_densities[q]
         * or @p material_model_outputs.densities[q], multiplied by the gravity vector
         * and dotted with the normal.
         * If the solid density is used, fluid is only fluxing in or out due to differences in
         * dynamic pressure, if the fluid density is used, melt fluxs in with the same velocity
         * as influxing solid material.
         *
         * @param boundary_indicator The boundary indicator of the part of the
         * boundary of the domain on which the point is located at which we
         * are requesting the fluid pressure gradients.
         * @param material_model_inputs The material property inputs.
         * @param material_model_outputs The material property outputs.
         * @return A vector of heatflux values at the evaluation points.
         */
        virtual
        std::vector<Tensor<1,dim> >
        heat_flux (const types::boundary_id boundary_indicator,
          const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
          const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
          const std::vector<Tensor<1,dim> > &normal_vectors) const = 0;

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
     * @ingroup BoundaryHeatFlux
     */
    template <int dim>
    void
    register_boundary_heat_flux (const std::string &name,
                                      const std::string &description,
                                      void (*declare_parameters_function) (ParameterHandler &),
                                      Interface<dim> *(*factory_function) ());

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
    Interface<dim> *
    create_boundary_heat_flux (ParameterHandler &prm);


    /**
     * Declare the runtime parameters of the registered fluid pressure boundary
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
     * for a fluid pressure boundary model, register it with the functions that
     * can declare their parameters and create these objects.
     *
     * @ingroup BoundaryHeatFlux
     */
#define ASPECT_REGISTER_BOUNDARY_HEAT_FLUX_MODEL(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_BOUNDARY_HEAT_FLUX_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::BoundaryHeatFlux::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::BoundaryHeatFlux::register_boundary_heat_flux<2>, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::BoundaryHeatFlux::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::BoundaryHeatFlux::register_boundary_heat_flux<3>, \
                                name, description); \
  }
  }
}


#endif
