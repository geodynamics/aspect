/*
  Copyright (C) 2013 - 2024 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/postprocess/visualization/dev_edot_ii.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      DevEdotii<dim>::
      DevEdotii ()
        :
        DataPostprocessorScalar<dim> ("dev_edot_ii",
                                      update_values | update_gradients | update_quadrature_points),
        Interface<dim>("1/s")
      {}



      template <int dim>
      void
      DevEdotii<dim>::
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

        MaterialModel::MaterialModelInputs<dim> in(input_data,
                                                    this->introspection());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
                                                      this->n_compositional_fields());

        this->get_material_model().create_additional_named_outputs(out);
        this->get_material_model().evaluate(in, out);        
        
        // This is only used for prescribed_dilation case
        MaterialModel::PrescribedPlasticDilation<dim>
        *prescribed_dilation = (this->get_parameters().enable_prescribed_dilation)
                              ? out.template get_additional_output<MaterialModel::PrescribedPlasticDilation<dim> >()
                              : nullptr;

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            SymmetricTensor<2,dim> strain_rate = in.strain_rate[q];

            // Note: Since the four cells share one interaction point (vertex), 
            // at the dike boundaries we only remove injection effects from points
            // where the injection rate is specified. To locate these prescribed
            // points, we artificially find points whose values greater than 0.9
            // times the prescribed values.
            if (prescribed_dilation != nullptr 
                && prescribed_dilation->dilation[q] != 0.0
                && this->get_timestep_number() > 0
                && std::fabs(strain_rate[0][0]) > 0.9 * prescribed_dilation->dilation[q])
             strain_rate[0][0] -= prescribed_dilation->dilation[q];

            computed_quantities[q](0) = std::sqrt(std::max(-second_invariant(deviator(strain_rate)), 0.));
          }

        // average the values if requested
        const auto &viz = this->get_postprocess_manager().template get_matching_postprocessor<Postprocess::Visualization<dim>>();
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(DevEdotii,
                                                  "deviatoric edot_ii",
                                                  "A visualization output object that generates output "
                                                  "for the norm of the deviatoric strain rate, i.e., for the quantity "
                                                  "$\\sqrt{\\varepsilon(\\mathbf u):\\varepsilon(\\mathbf u)}$ "
                                                  "in the incompressible case and "
                                                  "$\\sqrt{[\\varepsilon(\\mathbf u)-\\tfrac 13(\\textrm{tr}\\;\\varepsilon(\\mathbf u))\\mathbf I]:"
                                                  "[\\varepsilon(\\mathbf u)-\\tfrac 13(\\textrm{tr}\\;\\varepsilon(\\mathbf u))\\mathbf I]}$ "
                                                  "in the compressible case."
                                                  "\n\n"
                                                  "Physical units: \\si{\\per\\second}.")
    }
  }
}
