/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/dev_stress_tensor.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      DevStressTensor<dim>::
      DevStressTensor ()
        :
        DataPostprocessorTensor<dim> ("dev_stress_tensor",
                                      update_values | update_gradients | update_quadrature_points),
        Interface<dim>("Pa")
      {}



      template <int dim>
      void
      DevStressTensor<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert ((computed_quantities[0].size() == Tensor<2,dim>::n_independent_components),
                ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,   ExcInternalError());
        Assert (input_data.solution_gradients[0].size() == this->introspection().n_components,  ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> in(input_data,
                                                   this->introspection());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
                                                     this->n_compositional_fields());

        // We need to compute the viscosity
        in.requested_properties = MaterialModel::MaterialProperties::viscosity;

        this->get_material_model().create_additional_named_outputs(out);
        this->get_material_model().evaluate(in, out);      
        
        // This is only used for prescribed_dilation case
        MaterialModel::PrescribedPlasticDilation<dim>
        *prescribed_dilation = (this->get_parameters().enable_prescribed_dilation)
                              ? out.template get_additional_output<MaterialModel::PrescribedPlasticDilation<dim> >()
                              : nullptr;
        
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            // Compressive stress is negative by the sign convention
            // used by the engineering community, and as input and used
            // internally by ASPECT.
            // Here, we change the sign of the stress to match the
            // sign convention used by the geoscience community.
            SymmetricTensor<2,dim> stress = in.pressure[q] * unit_symmetric_tensor<dim>();
            
            // Initialization of the deviatoric stress
            SymmetricTensor<2,dim> deviatoric_sress = stress;

            // If elasticity is enabled, the deviatoric stress is stored
            // in compositional fields, otherwise the deviatoric stress
            // can be obtained from the viscosity and strain rate.
            if (this->get_parameters().enable_elasticity)
              {
                stress[0][0] -= in.composition[q][this->introspection().compositional_index_for_name("ve_stress_xx")];
                stress[1][1] -= in.composition[q][this->introspection().compositional_index_for_name("ve_stress_yy")];
                stress[0][1] -= in.composition[q][this->introspection().compositional_index_for_name("ve_stress_xy")];

                if (dim == 3)
                  {
                    stress[2][2] -= in.composition[q][this->introspection().compositional_index_for_name("ve_stress_zz")];
                    stress[0][2] -= in.composition[q][this->introspection().compositional_index_for_name("ve_stress_xz")];
                    stress[1][2] -= in.composition[q][this->introspection().compositional_index_for_name("ve_stress_yz")];
                  }
                
                deviatoric_sress = deviator(stress);
              }
            else
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

                const SymmetricTensor<2,dim> deviatoric_strain_rate = deviator(strain_rate);

                const double eta = out.viscosities[q];
                deviatoric_sress -= 2. * eta * deviatoric_strain_rate;
              }

            for (unsigned int d=0; d<dim; ++d)
              for (unsigned int e=0; e<dim; ++e)
                computed_quantities[q][Tensor<2,dim>::component_to_unrolled_index(TableIndices<2>(d,e))]
                  = deviatoric_sress[d][e];
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(DevStressTensor,
                                                  "deviatoric stress tensor",
                                                  "A visualization output object that generates output "
                                                  "for the 3 (in 2d) or 6 (in 3d) components of the deviatroic stress "
                                                  "tensor, i.e., for the components of the tensor "
                                                  "Physical units: \\si{\\pascal}.")
    }
  }
}
