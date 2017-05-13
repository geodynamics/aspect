#ifndef ASPECT_INCLUSION_COMPOSITIONAL_FIELDS_H
#define ASPECT_INCLUSION_COMPOSITIONAL_FIELDS_H

#include "../inclusion.h"



namespace aspect
{
  namespace InclusionBenchmark
  {
    using namespace dealii;

    template<int dim>
    class InclusionCompositionalMaterial : public InclusionMaterial<dim>
    {
      public:
        virtual double density(const double temperature,
                               const double pressure,
                               const std::vector<double> &compositional_fields,
                               const Point<dim> &position) const
        {
          return compositional_fields[0];
        }

        virtual double viscosity(const double temperature,
                                 const double pressure,
                                 const std::vector<double> &compositional_fields,
                                 const SymmetricTensor<2, dim> &strain_rate,
                                 const Point<dim> &position) const
        {
          return compositional_fields[1];
        }
    };
  }
}
#endif
