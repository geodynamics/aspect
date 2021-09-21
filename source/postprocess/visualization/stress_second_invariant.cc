/*
  Copyright (C) 2011 - 2021 by the authors of the ASPECT code.

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
                                      update_values | update_gradients | update_quadrature_points)
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

        MaterialModel::MaterialModelInputs<dim> in(input_data,
                                                   this->introspection());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
                                                     this->n_compositional_fields());

        this->get_material_model().create_additional_named_outputs(out);

        // Get the viscosity
        this->get_material_model().evaluate(in, out);

        for (unsigned int q = 0; q < n_quadrature_points; ++q)
          {
            const SymmetricTensor<2, dim> strain_rate = in.strain_rate[q];
            const SymmetricTensor<2, dim> deviatoric_strain_rate = (this->get_material_model().is_compressible()
                                                                    ? strain_rate - 1. / 3 * trace(strain_rate) * unit_symmetric_tensor<dim>()
                                                                    : strain_rate);

            const double eta = out.viscosities[q];

            // Compressive stress is positive in geoscience applications
            // TODO what sign to use here? In the material model we don't add
            // the minus sign for the comparison with the yield stress.
            SymmetricTensor<2, dim> stress = -2. * eta * deviatoric_strain_rate +
                                             in.pressure[q] * unit_symmetric_tensor<dim>();

            // Add elastic stresses if existent
            if (this->get_parameters().enable_elasticity == true)
              {
                stress[0][0] += in.composition[q][this->introspection().compositional_index_for_name("ve_stress_xx")];
                stress[1][1] += in.composition[q][this->introspection().compositional_index_for_name("ve_stress_yy")];
                stress[0][1] += in.composition[q][this->introspection().compositional_index_for_name("ve_stress_xy")];

                if (dim == 3)
                  {
                    stress[2][2] += in.composition[q][this->introspection().compositional_index_for_name("ve_stress_zz")];
                    stress[0][2] += in.composition[q][this->introspection().compositional_index_for_name("ve_stress_xz")];
                    stress[1][2] += in.composition[q][this->introspection().compositional_index_for_name("ve_stress_yz")];
                  }
              }

            // Compute the deviatoric stress after adding the elastic stresses that might have added
            // diagonal components
            const SymmetricTensor<2, dim> deviatoric_stress = deviator(stress);

            // Compute the second moment invariant of the deviatoric stress
            const double stress_invariant = std::sqrt(std::fabs(second_invariant(deviatoric_stress)));

            // Output the second stress invariant
            computed_quantities[q](0) = stress_invariant; 
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(StressSecondInvariant,
                                                  "stress second invariant",
                                                  "A visualization output object that generates output "
                                                  "the second moment invariant "
                                                  "of the deviatoric stress tensor.")
    }
  }
}
