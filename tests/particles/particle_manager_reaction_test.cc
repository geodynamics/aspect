// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

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
