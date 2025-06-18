/*
  Copyright (C) 2025 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/entropy_average.h>
#include <aspect/adiabatic_conditions/interface.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      EntropyAverage<dim>::
      EntropyAverage ()
        :
        DataPostprocessorScalar<dim> ("entropy_average",
                                      update_values | update_quadrature_points),
        Interface<dim>("J/kg/K")
      {}



      template <int dim>
      void
      EntropyAverage<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,           ExcInternalError());

        const std::vector<unsigned int> &chemical_composition_idx = this->introspection().chemical_composition_field_indices();

        const std::vector<unsigned int> &entropy_field_idx = this->introspection().get_indices_for_fields_of_type(CompositionalFieldDescription::entropy);

        const auto &composition_component_indices = this->introspection().component_indices.compositional_fields;

        Assert(chemical_composition_idx.size() == entropy_field_idx.size() - 1, ExcInternalError());

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            double sum_chemical_composition = 0.0;

            for (const auto chemical_composition_index: chemical_composition_idx)
              {
                sum_chemical_composition += input_data.solution_values[q][composition_component_indices[chemical_composition_index]];
              }

            sum_chemical_composition = std::min(sum_chemical_composition, 1.0);

            const double background_chemical_composition = 1.0 - sum_chemical_composition;

            double entropy = background_chemical_composition * input_data.solution_values[q][composition_component_indices[entropy_field_idx[0]]];

            for (unsigned int i=0; i<chemical_composition_idx.size(); ++i)
              {
                const unsigned int composition_component_index = composition_component_indices[chemical_composition_idx[i]];
                const unsigned int entropy_component_index = composition_component_indices[entropy_field_idx[i+1]];
                entropy += input_data.solution_values[q][composition_component_index] * input_data.solution_values[q][entropy_component_index];
              }

            computed_quantities[q](0) = entropy;
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(EntropyAverage,
                                                  "entropy average",
                                                  "A visualization output object that generates output "
                                                  "for the averaged entropy of multiple components."
                                                  "\n\n"
                                                  "Physical units: \\si{\\J/kg/K}.")
    }
  }
}
