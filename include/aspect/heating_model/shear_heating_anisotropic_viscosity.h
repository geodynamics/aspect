// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

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
