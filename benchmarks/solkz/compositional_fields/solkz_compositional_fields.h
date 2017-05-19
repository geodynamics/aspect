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
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const
        {
          SolKzMaterial<dim>::evaluate(in, out);

          for (unsigned int i=0; i < in.position.size(); ++i)
            {
              out.densities[i] = in.composition[i][0];
              out.viscosities[i] = in.composition[i][1];
            }
        }
    };
  }
}
#endif
