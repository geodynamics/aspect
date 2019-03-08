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
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const
        {
          for (unsigned int i=0; i < in.position.size(); ++i)
            {
              out.viscosities[i] = in.composition[i][0];
              out.densities[i] = 0;
              out.compressibilities[i] = 0;
              out.specific_heat[i] = 0;
              out.thermal_expansion_coefficients[i] = 0;
              out.thermal_conductivities[i] = 0.0;
            }
        }
    };
  }
}
#endif
