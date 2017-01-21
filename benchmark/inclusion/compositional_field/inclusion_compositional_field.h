#ifndef _aspect_inclusion_compositional_h
#define _aspect_inclusion_compositional_h

#include <../benchmark/inclusion/inclusion.h>

namespace aspect
{
  namespace InclusionBenchmark
  {
    using namespace dealii;

    template<int dim>
    class InclusionCompositionalMaterial : public InclusionMaterial<dim>
    {
      public:
        virtual double viscosity(const double temperature,
                                 const double pressure,
                                 const std::vector<double> &compositional_fields,
                                 const SymmetricTensor<2, dim> &strain_rate,
                                 const Point <dim> &position) const;
    };
  }
}
#endif
