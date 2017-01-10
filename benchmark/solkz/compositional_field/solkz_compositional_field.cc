#include "solkz_compositional_field.h"
#include <solkz.cc>

namespace aspect
{
  namespace InclusionBenchmark
  {
    using namespace dealii;

    template <int dim>
    double
    SolKzCompositionalMaterial<dim>::
    viscosity (const double,
               const double,
               const std::vector<double> &composition,
               const SymmetricTensor<2,dim> &,
               const Point<dim> &p) const
    {
      // defined as given in the Duretz et al. paper:
      // std::exp(2*0.5 * std::log(1e6)*p[1]);
      return composition[1];
    }

    template <int dim>
    double
    SolKzCompositionalMaterial<dim>::
    density (const double,
             const double,
             const std::vector<double> &composition,
             const Point<dim> &p) const
    {
      // defined as given in the paper
      // -std::sin(2*p[1])*std::cos(3*numbers::PI*p[0]);
      return composition[0];
    }
  }
}

// explicit instantiation
namespace aspect
{
  namespace InclusionBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(SolKzCompositionalMaterial,
                                   "SolKzCompositionalMaterial",
                                   "A material model that corresponds to the 'SolKz' benchmark "
                                   "using compositional fields as defined in Duretz et al., G-Cubed, 2011.")
  }
}
