/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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

#include <aspect/geometry_model/interface.h>
#include <aspect/material_model/composition_reaction.h>

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    class ParticleReactionTest : public MaterialModel::CompositionReaction<dim>
    {
      public:
        virtual void evaluate(const MaterialModelInputs<dim> &in,
                              MaterialModelOutputs<dim> &out) const
        {

          // Fill density, viscosity, conductivity, etc.
          CompositionReaction<dim>::evaluate(in, out);
          for (unsigned int q = 0; q < in.n_evaluation_points(); ++q)
            {
              // Remove the reaction produced by the base model.
              std::fill(out.reaction_terms[q].begin(),
                        out.reaction_terms[q].end(),
                        0.0);

              if (this->get_timestep_number() == 1)
                {
                  const double reaction_change = in.composition[q][0];
                  out.reaction_terms[q][0] -= reaction_change;
                  out.reaction_terms[q][1] += reaction_change;
                }
            }
        }
    };
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(ParticleReactionTest,
                                   "particle reaction test",
                                   "This is a simple material model to test "
                                   "'general composition reaction' particle "
                                   "property. Here it performs a reaction "
                                   "between two compositional fields at a "
                                   "specific timestep.");
  }
}
