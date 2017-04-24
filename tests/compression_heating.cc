#include <aspect/boundary_fluid_pressure/interface.h>
#include <aspect/material_model/melt_global.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/melt.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    template <int dim>
    class Compression : public MaterialModel::MeltGlobal<dim>
    {
      public:

        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;
    };


    template <int dim>
    void
    Compression<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      MeltGlobal<dim>::evaluate(in, out);
      const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
      for (unsigned int i=0; i < in.position.size(); ++i)
        {
          out.reaction_terms[i][porosity_idx] = 0.2;
        }
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Compression,
                                   "melt global compression",
                                   "A simple material model that is like the "
                                   "'melt simple' model, but has a constant reaction term.")
  }
}

