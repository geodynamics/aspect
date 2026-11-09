// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2017-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------



#ifndef _aspect_heating_model_compositional_heating_h
#define _aspect_heating_model_compositional_heating_h

#include <aspect/simulator_access.h>
#include <aspect/heating_model/interface.h>

namespace aspect
{
  namespace HeatingModel
  {
    /**
     * A class that implements a heating model where each compositional field
     * is assigned a user-defined internal heating value.
     *
     * @ingroup HeatingModels
     */
    template <int dim>
    class CompositionalHeating : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the magnitude of heat production arising from values assigned
         * to each compositional field.
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
         * Magnitude of heat production for each compositional field
         *
         * This variable is read from the parameter file through a parameter called 'Compositional heating values'.
         */
        std::vector<double> heating_values;

        /**
         * A vector with as many entries as compositional fields.
         * If the entry for a particular field is true the field is considered
         * during the averaging of heat production rates, if not the field is
         * ignored. This is useful if some compositional fields
         * are used to track properties like finite strain that should not
         * contribute to heat production.
         */
        std::vector<bool> fields_used_in_heat_production_averaging;

        /**
         * Similar to fields_used_in_heat_production_averaging, except it
         * determines whether to include the background field in the heat
         * production computation.
         */
        bool use_background_field_for_heat_production_averaging;
    };
  }
}


#endif
