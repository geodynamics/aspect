/*
  Copyright (C) 2024 by the authors of the ASPECT code.
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

#ifndef _aspect_heating_model_prescribed_dike_injection_h
#define _aspect_heating_model_prescribed_dike_injection_h

#include <aspect/heating_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  /* Head file for injection latent heat term*/
  namespace HeatingModel
  {
    using namespace dealii;

    /**
     * A class that implements the latent heat released during crystallization
     * of the melt lens and heating by melt injection into the model. It takes
     * the amount of material added on the right-hand side of the Stokes equations
     * and adds the corresponding heating term to the energy equation (considering
     * the latent heat of crystallization and the different temperature of the
     * injected melt).
     *
     * @ingroup HeatingModels
     */
    template <int dim>
    class LatentHeatDikeInjection : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Compute the heating model outputs for this class.
         */
        void
        evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                  const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                  HeatingModel::HeatingModelOutputs &heating_model_outputs) const override;
        void
        create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &material_model_outputs) const override;
        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */

        /**
         * Declare the parameters this class takes through input files.
         */
        static void declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void parse_parameters (ParameterHandler &prm) override;
        /**
         * @}
         */

      private:
        /**
         * Properties of injected material.
         */
        double latent_heat_of_crystallization;
        double temperature_of_injected_material;

        /**
         * Amount of new injected material from the dike
         */
        double dike_material_injection_fraction;
    };
  }
}

#endif