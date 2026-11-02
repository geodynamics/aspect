// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2015-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


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
         * Specify which material model outputs the heating model requires
         * for computing the heating terms.
         */
        MaterialModel::MaterialProperties::Property
        get_required_properties () const override;

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
        /**
         * Entropy change upon melting.
         *
         * This variable is read from the parameter file through a parameter called 'Melting entropy change'.
         */
        double melting_entropy_change;

        /**
         * This variable is read from the parameter file through a parameter called 'Retrieve entropy change from material model'.
         */
        bool   retrieve_entropy_change_from_material_model;
    };
  }
}


#endif
