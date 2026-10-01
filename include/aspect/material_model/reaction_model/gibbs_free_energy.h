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

#ifndef _aspect_material_model_reaction_model_gibbs_free_energy_h
#define _aspect_material_model_reaction_model_gibbs_free_energy_h

#include <aspect/material_model/reaction_model/gibbs_free_energy/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      template <int dim>
      struct GibbsFreeEnergyStep
      {
        std::shared_ptr<GibbsFreeEnergyModel::Interface<dim>> model;
        unsigned int local_reaction_index = numbers::invalid_unsigned_int;
      };

      template <int dim>
      class GibbsFreeEnergy
      {
        public:
          /**
           * Return the Gibbs free-energy difference G_B-G_A for a
           * reaction connecting adjacent phases A and B. The labels A and B
           * define the order of the difference and do not preassign the
           * instantaneous reactant and product.
           *
           * If G_B-G_A is negative, A -> B is favored, phase A is consumed,
           * and phase B is produced. If G_B-G_A is positive, B -> A is
           * favored, phase B is consumed, and phase A is produced. A zero
           * difference represents equilibrium. Phase availability and its
           * effect on reaction rates are handled by the kinetics model.
           */
          double delta_gibbs_free_energy(const double temperature,
                                         const double pressure,
                                         const GibbsFreeEnergyInputs &reference_state,
                                         const unsigned int reaction_index) const;

          unsigned int n_reactions() const;

          static void declare_parameters(ParameterHandler &prm);
          void parse_parameters(ParameterHandler &prm, const unsigned int n_reactions);

        private:
          std::vector<GibbsFreeEnergyStep<dim>> reactions;
      };
    }
  }
}

#endif
