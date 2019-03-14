#include <aspect/material_model/simple.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    template <int dim>
    class PrescribedFieldMaterial : public MaterialModel::Simple<dim>
    {
      public:

        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;

        virtual
        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const;
    };

  }
}

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    void
    PrescribedFieldMaterial<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      Simple<dim>::evaluate(in, out);

      // set up variable to interpolate prescribed field outputs onto compositional fields
      PrescribedFieldOutputs<dim> *prescribed_field_out = out.template get_additional_output<PrescribedFieldOutputs<dim> >();

      if (prescribed_field_out != NULL)
        for (unsigned int i=0; i < in.position.size(); ++i)
          prescribed_field_out->prescribed_field_outputs[i][1] = out.densities[i];
    }


    template <int dim>
    void
    PrescribedFieldMaterial<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output<PrescribedFieldOutputs<dim> >() == NULL)
        {
          const unsigned int n_points = out.viscosities.size();
          out.additional_outputs.push_back(
            std::make_shared<MaterialModel::PrescribedFieldOutputs<dim>> (n_points, this->n_compositional_fields()));
        }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(PrescribedFieldMaterial,
                                   "prescribed field material",
                                   "A simple material model that is like the "
                                   "'Simple' model, but creates prescribed field outputs.")
  }
}
