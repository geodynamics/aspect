/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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



#include <aspect/postprocess/visualization/stress_second_invariant.h>

#include <aspect/material_model/rheology/elasticity.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      StressSecondInvariant<dim>::
      StressSecondInvariant ()
        :
        DataPostprocessorScalar<dim> ("stress_second_invariant",
                                      update_values | update_gradients | update_quadrature_points),
        Interface<dim>("Pa")
      {}



      template <int dim>
      void
      StressSecondInvariant<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert(computed_quantities.size() == n_quadrature_points, ExcInternalError());
        Assert(computed_quantities[0].size() == 1, ExcInternalError());
        Assert(input_data.solution_values[0].size() == this->introspection().n_components, ExcInternalError());
        Assert(input_data.solution_gradients[0].size() == this->introspection().n_components, ExcInternalError());

        // Create the material model inputs and outputs to
        // retrieve the current viscosity.
        MaterialModel::MaterialModelInputs<dim> in(input_data,
                                                   this->introspection());

        in.requested_properties = MaterialModel::MaterialProperties::viscosity | MaterialModel::MaterialProperties::additional_outputs;

        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
                                                     this->n_compositional_fields());

        this->get_material_model().create_additional_named_outputs(out);

        this->get_material_model().evaluate(in, out);

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
          {
            // Compressive stress is negative by the sign convention
            // used by the engineering community, and as input and used
            // internally by ASPECT.
            // Here, we change the sign of the stress to match the
            // sign convention used by the geoscience community.
            // TODO Instead of first assembling the total stress
            // and then taking the deviator, we could only assemble
            // the deviatoric stress.
            SymmetricTensor<2,dim> stress = in.pressure[q] * unit_symmetric_tensor<dim>();

            // If elasticity is enabled, the deviatoric stress is stored
            // in compositional fields, otherwise the deviatoric stress
            // can be obtained from the viscosity and strain rate.
            if (this->get_parameters().enable_elasticity == true)
              {
                // Get the total deviatoric stress from the material model.
                const std::shared_ptr<const MaterialModel::ElasticAdditionalOutputs<dim>> elastic_additional_out
                  = out.template get_additional_output_object<MaterialModel::ElasticAdditionalOutputs<dim>>();

                Assert(elastic_additional_out != nullptr, ExcMessage("Elastic Additional Outputs are needed for the 'principal stress' postprocessor, but they have not been created."));

                stress -= elastic_additional_out->deviatoric_stress[q];
              }
            else
              {
                const SymmetricTensor<2,dim> strain_rate = in.strain_rate[q];
                const SymmetricTensor<2,dim> deviatoric_strain_rate
                  = (this->get_material_model().is_compressible()
                     ?
                     strain_rate - 1./3 * trace(strain_rate) * unit_symmetric_tensor<dim>()
                     :
                     strain_rate);

                const double eta = out.viscosities[q];
                stress -= 2. * eta * deviatoric_strain_rate;
              }

            // Compute the deviatoric stress tensor after elastic stresses were added.
            const SymmetricTensor<2, dim> deviatoric_stress = deviator(stress);


            // Compute the second moment invariant of the deviatoric stress
            // in the same way as the second moment invariant of the deviatoric
            // strain rate is computed in the viscoplastic material model.
            // TODO check that this is valid for the compressible case.
            const double stress_invariant = std::sqrt(std::max(-second_invariant(deviatoric_stress), 0.));

            computed_quantities[q](0) = stress_invariant;
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(StressSecondInvariant,
                                                  "stress second invariant",
                                                  "A visualization output object that outputs "
                                                  "the second moment invariant of the deviatoric stress tensor."
                                                  "\n\n"
                                                  "Physical units: $\\text{Pa}$.")
    }
  }
}
