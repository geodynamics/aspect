/*
  Copyright (C) 2022 - 2024 by the authors of the ASPECT code.

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

#include <aspect/material_model/composition_reaction.h>
#include <aspect/geometry_model/interface.h>

/**
 * This material model assumes three compositional fields
 * where the reaction rate of the first and the third one
 * depend on the second one. Thus, an iterated scheme is
 * required to compute the correct solution.
 */

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    class IteratedReaction : public MaterialModel::CompositionReaction<dim>
    {
      public:
        virtual void evaluate(const MaterialModelInputs<dim> &in,
                              MaterialModelOutputs<dim> &out) const
        {
          this->CompositionReaction<dim>::evaluate(in, out);
          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
            {
              const double depth = this->get_geometry_model().depth(in.position[i]);
              for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                {
                  Assert(in.composition[i].size() > 1,
                         ExcMessage ("Material model iterated reaction can only be used with "
                                     "at least two compositial fields."));

                  double delta_C = 0.0;
                  switch (c)
                    {
                      case 0:
                        delta_C = in.composition[i][1];
                        break;
                      case 1:
                        delta_C = 1.0;
                        break;
                      case 2:
                        delta_C = in.composition[i][1];
                        break;
                    }
                  out.reaction_terms[i][c] = delta_C;
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
    ASPECT_REGISTER_MATERIAL_MODEL(IteratedReaction,
                                   "iterated reaction",
                                   "A simple material model that is like the "
                                   "'composition reaction' model, but requires an "
                                   "iterated Advection and Stokes scheme to converge to the correct "
                                   "solution.")
  }
}
