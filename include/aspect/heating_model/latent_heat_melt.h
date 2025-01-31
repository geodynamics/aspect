/*
  Copyright (C) 2015 - 2022 by the authors of the ASPECT code.

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


#ifndef _aspect_heating_model_latent_heat_melt_h
#define _aspect_heating_model_latent_heat_melt_h

#include <aspect/heating_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace HeatingModel
  {
    /**
     * A class that implements a standard formulation of latent heat
     * of melting. This assumes that there is a compositional field
     * called porosity, and it uses the reaction term of this field
     * (the fraction of material that melted in the current time step)
     * multiplied by a constant entropy change for melting 100%
     * of the material as source term of the heating model.
     * The left-hand side term is zero.
     *
     * @ingroup HeatingModels
     */
    template <int dim>
    class LatentHeatMelt : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Compute the heating model outputs for this class.
         */
        void
        evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                  const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                  HeatingModel::HeatingModelOutputs &heating_model_outputs) const override;

        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */
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

        /**
         * @}
         */

        virtual
        void
        create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const override;

      private:
        // entropy change upon melting
        double melting_entropy_change;
        bool   retrieve_entropy_change_from_material_model;
    };
  }
}


#endif
