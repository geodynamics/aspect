#include "solcx_compositional_field.h"
#include <../benchmark/solcx/solcx.cc>

namespace aspect
{
  namespace InclusionBenchmark
  {
    using namespace dealii;

    template <int dim>
    double
    SolCxCompositionalMaterial<dim>::
    viscosity (const double,
               const double,
               const std::vector<double> &composition,
               const SymmetricTensor<2,dim> &,
               const Point<dim> &p) const
    {
      // defined as given in the Duretz et al. paper: (p[0] < 0.5 ? 1 : eta_B)
      // (p[0] < 0.5 ? 1 : eta_B)
      return composition[1];
    }

    template <int dim>
    double
    SolCxCompositionalMaterial<dim>::
    density (const double,
             const double,
             const std::vector<double> &composition,
             const Point<dim> &p) const
    {
      // defined as given in the paper, plus the constant background density:
      // background_density - std::sin(numbers::PI*p[1])*std::cos(numbers::PI*p[0])
      return composition[0];
    }
  }
}



// explicit instantiations
namespace aspect
{
  namespace InclusionBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(SolCxCompositionalMaterial,
                                   "SolCxCompositionalMaterial",
                                   "A material model that corresponds to the 'SolCx' benchmark "
                                   "defined in Duretz et al., G-Cubed, 2011.")
  }
}
