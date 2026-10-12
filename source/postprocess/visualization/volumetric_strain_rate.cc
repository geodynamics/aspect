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


#include <aspect/postprocess/visualization/volumetric_strain_rate.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      VolumetricStrainRate<dim>::
      VolumetricStrainRate ()
        :
        DataPostprocessorScalar<dim> ("volumetric_strain_rate",
                                      update_gradients | update_quadrature_points),
        Interface<dim>("1/s")
      {}



      template <int dim>
      void
      VolumetricStrainRate<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,           ExcInternalError());
        Assert (input_data.solution_gradients[0].size() == this->introspection().n_components,          ExcInternalError());

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            // extract the velocity gradients
            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = input_data.solution_gradients[q][d];

            const double div_u = trace(grad_u);
            computed_quantities[q](0) = div_u;
          }

        // average the values if requested
        const auto &viz = this->get_postprocess_manager().template get_matching_active_plugin<Postprocess::Visualization<dim>>();
        if (!viz.output_pointwise_stress_and_strain())
          average_quantities(computed_quantities);
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(VolumetricStrainRate,
                                                  "volumetric strain rate",
                                                  "A visualization output object that generates output "
                                                  "for the volumetric strain rate, i.e., for the quantity "
                                                  "$\\nabla\\cdot\\mathbf u = \\textrm{div}\\; \\mathbf u = "
                                                  "\\textrm{trace}\\; \\varepsilon(\\mathbf u)$. "
                                                  "This should be zero (in some average sense) in incompressible "
                                                  "convection models, but can be non-zero in compressible models and "
                                                  "models with melt transport."
                                                  "\n\n"
                                                  "Physical units: $\\frac{1}{\\text{s}}$.")
    }
  }
}
