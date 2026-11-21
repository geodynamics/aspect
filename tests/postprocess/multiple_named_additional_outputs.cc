// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2022-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include <aspect/material_model/simple.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
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
      const std::shared_ptr<PrescribedFieldOutputs<dim>> prescribed_field_out
        = out.template get_additional_output_object<PrescribedFieldOutputs<dim>>();

      if (prescribed_field_out != nullptr)
        for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
          {
            const double y = in.position[i](1);
            prescribed_field_out->prescribed_field_outputs[i][1] = std::exp(-y*y/2.0);
          }

      const std::shared_ptr<SeismicAdditionalOutputs<dim>> seismic_out
        = out.template get_additional_output_object<SeismicAdditionalOutputs<dim>>();

      if (seismic_out != nullptr)
        for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
          {
            const double y = in.position[i](1);
            seismic_out->vp[i] = std::exp(-y*y/3.0);
            seismic_out->vs[i] = -std::exp(-y*y/3.0);
          }
    }


    template <int dim>
    void
    PrescribedFieldMaterial<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template has_additional_output_object<PrescribedFieldOutputs<dim>>() == false)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::PrescribedFieldOutputs<dim>> (n_points,this->n_compositional_fields()));
        }

      if (out.template has_additional_output_object<SeismicAdditionalOutputs<dim>>() == false)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::SeismicAdditionalOutputs<dim>> (n_points));
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
