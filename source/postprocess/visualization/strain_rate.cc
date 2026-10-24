// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/postprocess/visualization/strain_rate.h>
#include <aspect/utilities.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      StrainRate<dim>::
      StrainRate ()
        :
        DataPostprocessorScalar<dim> ("strain_rate",
                                      update_gradients | update_quadrature_points),
        Interface<dim>("1/s")
      {}



      template <int dim>
      void
      StrainRate<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,
                ExcInternalError());
        Assert (input_data.solution_gradients[0].size() == this->introspection().n_components,
                ExcInternalError());

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            Tensor<2,dim> grad_u;
            for (unsigned int d=0; d<dim; ++d)
              grad_u[d] = input_data.solution_gradients[q][d];

            const SymmetricTensor<2,dim> strain_rate = symmetrize (grad_u);
            computed_quantities[q](0) = std::sqrt(std::max(-Utilities::Tensors::consistent_second_invariant_of_deviatoric_tensor(Utilities::Tensors::consistent_deviator(strain_rate)), 0.));
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(StrainRate,
                                                  "strain rate",
                                                  "A visualization output object that generates output "
                                                  "for the norm of the strain rate, i.e., for the quantity "
                                                  "$\\sqrt{\\varepsilon(\\mathbf u):\\varepsilon(\\mathbf u)}$ "
                                                  "in the incompressible case and "
                                                  "$\\sqrt{[\\varepsilon(\\mathbf u)-\\tfrac 13(\\textrm{tr}\\;\\varepsilon(\\mathbf u))\\mathbf I]:"
                                                  "[\\varepsilon(\\mathbf u)-\\tfrac 13(\\textrm{tr}\\;\\varepsilon(\\mathbf u))\\mathbf I]}$ "
                                                  "in the compressible case. It is also called effective deviatoric strain rate "
                                                  "in Glerum et al. (2018)."
                                                  "\n\n"
                                                  "Physical units: $\\frac{1}{\\text{s}}$.")
    }
  }
}
