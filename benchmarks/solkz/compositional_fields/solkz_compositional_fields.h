#ifndef ASPECT_SOLKZ_COMPOSITIONAL_FIELDS_H
#define ASPECT_SOLKZ_COMPOSITIONAL_FIELDS_H


#include "../solkz.h"



namespace aspect
{
  namespace InclusionBenchmark
  {
    using namespace dealii;

    template<int dim>
    class SolKzCompositionalMaterial : public SolKzMaterial<dim>
    {
      public:
        virtual double density(const double temperature,
                               const double pressure,
                               const std::vector<double> &composition,
                               const Point<dim> &position) const
        {
          return composition[0];
        }

        virtual double viscosity(const double temperature,
                                 const double pressure,
                                 const std::vector<double> &composition,
                                 const SymmetricTensor<2, dim> &strain_rate,
                                 const Point<dim> &position) const
        {
          return composition[1];
        }
    };
  }
}
#endif
