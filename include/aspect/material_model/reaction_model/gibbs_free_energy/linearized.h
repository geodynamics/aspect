/*
  Copyright (C) 2011 - 2026 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_reaction_model_gibbs_free_energy_linearized_h
#define _aspect_material_model_reaction_model_gibbs_free_energy_linearized_h

#include <aspect/material_model/reaction_model/gibbs_free_energy/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      namespace GibbsFreeEnergyModel
      {
        template <int dim>
        class Linearized : public Interface<dim>
        {
          public:
            double compute_delta_gibbs_free_energy(const double temperature,
                                                   const double pressure,
                                                   const GibbsFreeEnergyInputs &reference_state,
                                                   const unsigned int reaction_index) const override;

            static void declare_parameters(ParameterHandler &prm);
            void parse_parameters(ParameterHandler &prm, const unsigned int n_reactions) override;
        };
      }
    }
  }
}

#endif
