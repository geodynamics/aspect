#include "inclusion_compositional_field.h"
#include <inclusion.cc>

namespace aspect
{
  namespace InclusionBenchmark
  {
    using namespace dealii;

    template <int dim>
    double
    InclusionCompositionalMaterial<dim>::
    viscosity (const double,
               const double,
               const std::vector<double> &composition,
               const SymmetricTensor<2,dim> &,
               const Point<dim> &p) const
    {
      return composition[0];
    }
  }
}



// explicit instantiations
namespace aspect
{
  namespace InclusionBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(InclusionCompositionalMaterial,
                                   "InclusionCompositionalMaterial",
                                   "A material model that corresponds to the 'Inclusion' benchmark "
                                   "defined in Duretz et al., G-Cubed, 2011 using compositional fields.")
  }
}
