/*
  Copyright (C) 2025 - by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file doc/COPYING. If not see
  <http://www.gnu.org/licenses/>.
*/

#ifndef _aspect_material_model_additional_outputs_peng_robinson_fugacity76_h
#define _aspect_material_model_additional_outputs_peng_robinson_fugacity76_h


#include <aspect/material_model/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * Additional output field named "fugacity" for the pure fluid configured
     * for the Peng-Robinson equation of state.
     */
    template <int dim>
    class PengRobinsonFugacity : public NamedAdditionalMaterialOutputs<dim>
    {
      public:
        /**
         * Constructor. Allocate output storage for @p n_points.
         */
        PengRobinsonFugacity(const unsigned int n_points);

        /**
         * Return the requested named output.
         */
        std::vector<double>
        get_nth_output(const unsigned int idx) const override;

        /**
         * Peng-Robinson fugacity of a pure fluid in Pa for each compositional field
         */
        std::vector<double> fugacities;
    };
  }
}

#endif
