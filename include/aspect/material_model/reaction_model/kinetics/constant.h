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

#ifndef _aspect_material_model_reaction_model_kinetics_constant_h
#define _aspect_material_model_reaction_model_kinetics_constant_h

#include <aspect/material_model/reaction_model/kinetics/cahn1956_interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      /**
       * A trivial kinetics model that returns a user-specified constant
       * reaction rate. This model is intended for testing reaction chains.
       */
      template <int dim>
      class ConstantReactionRate : public Cahn1956Interface<dim>
      {
        public:
          double net_forward_reaction_rate(const double temperature,
                                           const double pressure,
                                           const double delta_forward_gibbs_energy,
                                           const double cumulative_forward_reaction_progress,
                                           const unsigned int reaction_index) const override;

          static void declare_parameters(ParameterHandler &prm);

          void parse_parameters(ParameterHandler &prm, const unsigned int n_reactions) override;

        private:
          /**
           * The user-provided constant reaction rates.
           */
          std::vector<double> reaction_rates;
      };
    }
  }
}

#endif
