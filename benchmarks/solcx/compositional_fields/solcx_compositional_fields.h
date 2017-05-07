#ifndef ASPECT_SOLCX_COMPOSITIONAL_FIELDS_H
#define ASPECT_SOLCX_COMPOSITIONAL_FIELDS_H

#include "../solcx.h"



namespace aspect
{
  namespace InclusionBenchmark
  {
    using namespace dealii;

    template<int dim>
    class SolCxCompositionalMaterial : public SolCxMaterial<dim>
    {
      public:
        double density(const double temperature,
                       const double pressure,
                       const std::vector<double> &composition,
                       const Point <dim> &position) const
        {
          return composition[0];
        }

        double viscosity(const double temperature,
                         const double pressure,
                         const std::vector<double> &composition,
                         const SymmetricTensor<2, dim> &strain_rate,
                         const Point <dim> &position) const
        {
          return composition[1];
        }
    };
  }
}
#endif
