/*
 Copyright (C) 2025 by the authors of the ASPECT code.

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

#ifndef _aspect_heating_model_shear_heating_anisotropic_viscosity_h
#define _aspect_heating_model_shear_heating_anisotropic_viscosity_h

#include <aspect/heating_model/interface.h>

#include <aspect/simulator_access.h>

namespace aspect
{
  namespace HeatingModel
  {
    /**
     * A class that implements a standard model for shear heating extended for an
     * anisotropic viscosity tensor. If the material model provides a stress-strain
     * director tensor, then the strain-rate is multiplied with this
     * tensor to compute the stress that is used when computing the shear heating.
     *
     * @ingroup HeatingModels
     */
    template <int dim>
    class ShearHeatingAnisotropicViscosity : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
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
         * Allow the heating model to attach additional material model outputs.
         */
        void
        create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &material_model_outputs) const override;
    };
  }
}

#endif
