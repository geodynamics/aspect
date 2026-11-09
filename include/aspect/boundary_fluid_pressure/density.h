// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2015-2021 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_boundary_fluid_pressure_density_h
#define _aspect_boundary_fluid_pressure_density_h

#include <aspect/boundary_fluid_pressure/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace BoundaryFluidPressure
  {
    /**
     * A class that implements simple fluid pressure boundary conditions
     * based on the fluid or solid density from the material model.
     *
     * The fluid pressure gradient can either be set to
     * $\rho_s \textbf{g}$ (solid density times gravity) or to
     * $\rho_f \textbf{g}$ (fluid density times gravity).
     *
     * @ingroup BoundaryFluidPressures
     */
    template <int dim>
    class Density : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * @copydoc Interface::fluid_pressure_gradient
         */
        void fluid_pressure_gradient (
          const types::boundary_id boundary_indicator,
          const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
          const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
          const std::vector<Tensor<1,dim>> &normal_vectors,
          std::vector<double> &fluid_pressure_gradient_outputs
        ) const override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * Identify which density to use to compute the fluid pressure
         * gradient at the model boundary.
         */
        struct DensityFormulation
        {
          enum Kind
          {
            solid_density,
            fluid_density,
            average_density
          };
        };

        typename DensityFormulation::Kind density_formulation;
    };
  }
}


#endif
