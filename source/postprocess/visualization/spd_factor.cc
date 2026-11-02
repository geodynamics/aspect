// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2017-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/postprocess/visualization/spd_factor.h>
#include <aspect/newton.h>
#include <aspect/utilities.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      SPD_Factor<dim>::
      SPD_Factor ()
        :
        DataPostprocessorScalar<dim> ("spd_factor",
                                      update_values | update_gradients | update_quadrature_points),
        Interface<dim>("")  // no physical units
      {}



      template <int dim>
      void
      SPD_Factor<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        Assert(Parameters<dim>::is_defect_correction(this->get_parameters().nonlinear_solver),
               ExcMessage("The SPD factor plugin can only be used with defect correction type Stokes or Newton Stokes "
                          "solvers."));

        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,                             ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                                            ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,    ExcInternalError());
        Assert (input_data.solution_gradients[0].size() == this->introspection().n_components, ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> in(input_data,
                                                   this->introspection());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
                                                     this->n_compositional_fields());
        in.requested_properties = MaterialModel::MaterialProperties::viscosity | MaterialModel::MaterialProperties::additional_outputs;

        out.additional_outputs.push_back(
          std::make_unique<MaterialModel::MaterialModelDerivatives<dim>> (n_quadrature_points));

        this->get_material_model().evaluate(in, out);

        const std::shared_ptr<const MaterialModel::MaterialModelDerivatives<dim>> derivatives
          = out.template get_additional_output_object<MaterialModel::MaterialModelDerivatives<dim>>();

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            computed_quantities[q](0) = Utilities::compute_spd_factor<dim>(out.viscosities[q],
                                                                           in.strain_rate[q],
                                                                           derivatives->viscosity_derivative_wrt_strain_rate[q],
                                                                           this->get_newton_handler().parameters.SPD_safety_factor);
          }
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(SPD_Factor,
                                                  "spd factor",
                                                  "A visualization output object that generates output "
                                                  "for the spd factor. The spd factor is a factor which "
                                                  "scales a part of the Jacobian used for the Newton solver "
                                                  "to make sure that the Jacobian remains positive definite."
                                                  "\n\n"
                                                  "Physical units: None.")
    }
  }
}
