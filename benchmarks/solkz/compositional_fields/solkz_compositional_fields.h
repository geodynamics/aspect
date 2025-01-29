#ifndef ASPECT_SOLKZ_COMPOSITIONAL_FIELDS_H
#define ASPECT_SOLKZ_COMPOSITIONAL_FIELDS_H


#include "../solkz.h"



namespace aspect
{
  namespace InclusionBenchmark
  {
    template <int dim>
    class SolKzCompositionalMaterial : public SolKzMaterial<dim>
    {
      public:
        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override
        {
          SolKzMaterial<dim>::evaluate(in, out);

          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
            {
              out.densities[i] = in.composition[i][0];
              out.viscosities[i] = in.composition[i][1];
            }
        }
    };
  }
}
#endif
