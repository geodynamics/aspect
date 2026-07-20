/*
  Copyright (C) 2015 - 2023 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/prescribed_solution.h>
#include <aspect/prescribed_solution/interface.h>
#include <aspect/utilities.h>
#include <aspect/material_model/interface.h>
#include <deal.II/numerics/data_out.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
      * A postprocessor that visualizes the analytical solution.
      */
      template <int dim>
      PrescribedSolutionPostprocessor<dim>::
      PrescribedSolutionPostprocessor()
        :
        DataPostprocessor<dim> ()
      {}



      template <int dim>
      std::vector<std::string>
      PrescribedSolutionPostprocessor<dim>::
      get_names () const
      {
        std::vector<std::string> solution_names;

        for (unsigned int component_index=0; component_index< this->introspection().n_components; ++component_index)
          {
            if (component_index == this->introspection().component_indices.temperature)
              {
                solution_names.push_back("prescribed_temperature_indicator");
                solution_names.push_back("prescribed_temperature_value");
              }
            else if (component_index == this->introspection().component_indices.pressure)
              {
                solution_names.push_back("prescribed_pressure_indicator");
                solution_names.push_back("prescribed_pressure_value");
              }
            else if ((component_index >= this->introspection().component_indices.velocities[0]) &&
                     (component_index <= this->introspection().component_indices.velocities[dim-1]))
              {
                solution_names.push_back("prescribed_velocity_indicator_" +
                                         std::to_string(component_index-this->introspection().component_indices.velocities[0]));
                solution_names.push_back("prescribed_velocity_value_" +
                                         std::to_string(component_index-this->introspection().component_indices.velocities[0]));
              }
            else if ((component_index >= this->introspection().component_indices.compositional_fields[0]) &&
                     (component_index <= this->introspection().n_components-1))
              {
                solution_names.push_back("prescribed_composition_indicator_" +
                                         std::to_string(component_index-this->introspection().component_indices.compositional_fields[0]));
                solution_names.push_back("prescribed_composition_value_" +
                                         std::to_string(component_index-this->introspection().component_indices.compositional_fields[0]));
              }
          }

        return solution_names;
      }



      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      PrescribedSolutionPostprocessor<dim>::
      get_data_component_interpretation () const
      {
        std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation;

        // indicators and values
        for (unsigned int component_index=0; component_index< this->introspection().n_components; ++component_index)
          {
            interpretation.push_back (DataComponentInterpretation::component_is_scalar);
            interpretation.push_back (DataComponentInterpretation::component_is_scalar);
          }

        return interpretation;
      }



      template <int dim>
      UpdateFlags
      PrescribedSolutionPostprocessor<dim>::
      get_needed_update_flags () const
      {
        return update_quadrature_points;
      }



      template <int dim>
      void
      PrescribedSolutionPostprocessor<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_evaluation_points = input_data.solution_values.size();

        Assert (this->get_parameters().include_melt_transport == false, ExcInternalError());
        Assert (computed_quantities.size() == n_evaluation_points, ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components, ExcInternalError());

        const aspect::PrescribedSolution::Manager<dim> &prescribed_solution_manager
          = this->get_prescribed_solution();
        const auto &plugin_objects = prescribed_solution_manager.get_active_plugins();

        // Variables to store constrained solution
        std::vector<unsigned int> component_indices(n_evaluation_points);
        std::vector<bool> component_is_constrained(n_evaluation_points);
        std::vector<double> constrained_component_value(n_evaluation_points);

        // determine prescribed solution components
        for (unsigned int component_index=0; component_index<this->introspection().n_components; ++component_index)
          {
            // assign initial 0 values
            for (unsigned int q=0; q<n_evaluation_points; ++q)
              {
                component_indices[q] = component_index;
                component_is_constrained[q] = false;
                constrained_component_value[q] = 0.0;

                computed_quantities[q][2*component_index] = 0.0;
                computed_quantities[q][2*component_index+1] = 0.0;
              }

            for (auto &it: plugin_objects)
              {
                typename DoFHandler<dim>::active_cell_iterator cell;

                it->constrain_solution(cell,
                                       input_data.evaluation_points,
                                       component_indices,
                                       component_is_constrained,
                                       constrained_component_value);

              }

            for (unsigned int q=0; q<n_evaluation_points; ++q)
              {
                computed_quantities[q][2*component_index] = std::max(computed_quantities[q][2*component_index], static_cast<double>(component_is_constrained[q]));

                if (component_is_constrained[q] == true)
                  computed_quantities[q][2*component_index+1] = constrained_component_value[q];
              }
          }
      }
    }
  }
}


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {

      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(PrescribedSolutionPostprocessor,
                                                  "prescribed solution",
                                                  "A visualization output object that outputs "
                                                  "whether the solution components are prescribed by "
                                                  "the prescribed solution plugin system and if so to "
                                                  "which value.")

    }
  }
}
