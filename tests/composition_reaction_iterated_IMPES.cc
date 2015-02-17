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
        /**
         * @name Physical parameters used in the basic equations
         * @{
         */
        virtual double reaction_term (const double temperature,
                                      const double pressure,
                                      const std::vector<double> &compositional_fields,
                                      const Point<dim> &position,
                                      const unsigned int compositional_variable) const;
    };

  }
}

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    double
    IteratedReaction<dim>::
    reaction_term (const double temperature,
                   const double pressure,
                   const std::vector<double> &compositional_fields,
                   const Point<dim> &position,
                   const unsigned int compositional_variable) const
    {
      Assert(compositional_fields.size() > 1, 
            ExcMessage ("Material model iterated reaction can only be used with "
                        "at least two compositial fields."));
      double delta_C = 0.0;
      switch (compositional_variable)
        {
          case 0:
            delta_C = compositional_fields[1];
            break;
          case 1:
            delta_C = 1;
            break;
          case 2:
            delta_C = compositional_fields[1];
            break;
          default:
            delta_C = 0.0;
            break;
        }
      return delta_C;
    }
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
