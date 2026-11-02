// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2015-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_heating_model_shear_heating_with_melt_h
#define _aspect_heating_model_shear_heating_with_melt_h

#include <aspect/heating_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace HeatingModel
  {
    /**
     * A class that implements a standard shear heating rate for dissipation
     * heating of partially molten material.
     *
     * Add the term
     *    $   \xi \left( \nabla \cdot \mathbf u_s \right)^2 $,
     * which corresponds to bulk (compression) heating,
     * and the term
     *    $   \frac{\eta_f \phi^2}{k_{\phi}} \left( \mathbf u_f - \mathbf u_s \right)^2 $,
     * which corresponds to heating due to melt segregation.
     * Both terms are heating source terms, and the left-hand side term is zero.
     *
     * Also see the Equations section in the manual.
     *
     * @ingroup HeatingModels
     */
    template <int dim>
    class ShearHeatingMelt : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
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
         * Specify which material model outputs the heating model requires
         * for computing the heating terms.
         */
        MaterialModel::MaterialProperties::Property
        get_required_properties () const override;

        /**
         * Allow the heating model to attach additional material model outputs.
         */
        void
        create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &) const override;

        /**
         * Allow the heating model to attach additional material model inputs it needs.
         */
        void
        create_additional_material_model_inputs(MaterialModel::MaterialModelInputs<dim> &) const override;
    };
  }
}


#endif
