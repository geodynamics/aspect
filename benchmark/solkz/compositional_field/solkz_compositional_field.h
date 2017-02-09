#ifndef _aspect_solkz_compositional_field_h
#define _aspect_solkz_compositional_field_h

#include <../benchmark/solkz/solkz.h>

namespace aspect
{
  namespace InclusionBenchmark
  {
    using namespace dealii;

    template<int dim>
    class SolKzCompositionalMaterial : public SolKzMaterial<dim>
    {
      public:
        virtual double viscosity(const double temperature,
                                 const double pressure,
                                 const std::vector<double> &compositional_fields,
                                 const SymmetricTensor<2, dim> &strain_rate,
                                 const Point <dim> &position) const;

        virtual double density(const double temperature,
                               const double pressure,
                               const std::vector<double> &compositional_fields,
                               const Point <dim> &position) const;
    };
  }
}
#endif
