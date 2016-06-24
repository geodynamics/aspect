#include <aspect/material_model/composition_reaction.h>

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
    using namespace dealii;

    template <int dim>
    class IteratedReaction : public MaterialModel::CompositionReaction<dim>
    {
      public:
        virtual void evaluate(const MaterialModelInputs<dim> &in,
                              MaterialModelOutputs<dim> &out) const
        {
          this->CompositionReaction<dim>::evaluate(in, out);
          for (unsigned int i=0; i < in.position.size(); ++i)
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
                                   "iterated IMPES scheme to converge to the correct "
                                   "solution.")
  }
}
