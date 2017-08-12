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
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const
        {
          for (unsigned int i=0; i < in.position.size(); ++i)
            {
              const Point<dim> &pos = in.position[i];
              out.densities[i] = in.composition[i][0];
              out.viscosities[i] = in.composition[i][1];
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
