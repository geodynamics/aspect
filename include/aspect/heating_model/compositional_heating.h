/*
  Copyright (C) 2017 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
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
    using namespace dealii;

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
         * Constructor.
         */
        CompositionalHeating ();

        /**
         * Return the magnitude of heat production arising from values assigned
         * to each compositional field.
         */
        virtual
        void
        evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                  const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                  HeatingModel::HeatingModelOutputs &heating_model_outputs) const;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:

        /**
         * From a list of compositional fields of length N, we come up with an
         * N+1 length list that which also includes the fraction of
         * ``background mantle''. This list should sum to one, and is
         * interpreted as volume fractions.  If the sum of the
         * compositional_fields is greater than one, we assume that there is
         * no background mantle (i.e., that field value is zero).  Otherwise,
         * the difference between the sum of the compositional fields and 1.0
         * is assumed to be the amount of background mantle.
         */
        std::vector<double> compute_volume_fractions(
          const std::vector<double> &compositional_fields) const;

        /**
         * Magnitude of heat production in each compositional field
         */
        std::vector<double> heating_values;

        /**
         * User can choose whether the heat production value associated with
         * each compositional field is used when the average value at each
         * point is calculated. This is useful if some compositional fields
         * are used to track properties like finite strain that should not
         * contribute to material properties like viscosity, heat production,
         * etc.
         */
        std::vector<int> field_used_in_heat_production_averaging;

    };
  }
}


#endif
