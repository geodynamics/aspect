/*
  Copyright (C) 2017 - 2019 by the authors of the ASPECT code.

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
