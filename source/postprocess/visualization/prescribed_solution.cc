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
#include <aspect/simulator.h>
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

        // for temperature indicators
        solution_names.push_back("prescribed_temperature_indicator");

        // for temperature values
        solution_names.push_back("prescribed_temperature");

        // for velocity indicators
        solution_names.push_back("prescribed_velocity_indicator");

        // for velocity values
        for (unsigned int d=0; d<dim; ++d)
          {
            solution_names.push_back("prescribed_velocity");
          }

        // for composition indicators
        if (this->n_compositional_fields() > 0)
          solution_names.push_back("prescribed_composition_indicator");

        // for composition values
        for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
          {
            solution_names.push_back("prescribed_composition_" + std::to_string(c));
          }

        return solution_names;
      }



      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      PrescribedSolutionPostprocessor<dim>::
      get_data_component_interpretation () const
      {
        std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation;

        // temperature indicators
        interpretation.push_back (DataComponentInterpretation::component_is_scalar);

        // temperature values
        interpretation.push_back (DataComponentInterpretation::component_is_scalar);

        // velocity indicators
        interpretation.push_back (DataComponentInterpretation::component_is_scalar);

        // velocity values
        for (unsigned int d=0; d<dim; ++d)
          interpretation.push_back (DataComponentInterpretation::component_is_part_of_vector);

        // composition indicators
        if (this->n_compositional_fields() > 0)
          interpretation.push_back (DataComponentInterpretation::component_is_scalar);

        // composition values
        for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
          {
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
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,           ExcInternalError());


        const aspect::PrescribedSolution::Manager<dim> &prescribed_solution_manager
          = Plugins::get_plugin_as_type<const aspect::PrescribedSolution::Manager<dim>>(this->get_prescribed_solution());

        const auto &plugin_objects = prescribed_solution_manager.get_active_plugins();

        const unsigned int n_dofs = input_data.evaluation_points.size();

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            // assign initial 0 values
            const unsigned int total_output = (this->n_compositional_fields() > 0)? 4+dim+this->n_compositional_fields(): 3+dim;
            for (unsigned int output_index = 0; output_index < total_output; ++output_index)
              computed_quantities[q][output_index] = 0.0;

            for (auto &it: plugin_objects)
              {
                typename DoFHandler<dim>::active_cell_iterator cell;

                // for temperature
                std::vector<unsigned int> temperature_component_indices(n_dofs, this->introspection().component_indices.temperature);
                std::vector<bool> temperature_should_be_constrained(n_dofs);
                std::vector<double> temperature_solution(n_dofs);
                it->constrain_solution(cell,
                                       input_data.evaluation_points,
                                       temperature_component_indices,
                                       temperature_should_be_constrained,
                                       temperature_solution);


                computed_quantities[q][0] = std::max(computed_quantities[q][0], double(temperature_should_be_constrained[q]));

                if (temperature_should_be_constrained[q] == true)
                  computed_quantities[q][1] = temperature_solution[q];

                // for velocity
                for (unsigned int d=0; d<dim; ++d)
                  {
                    // specify the component
                    std::vector<unsigned int> velocity_component_indices(n_dofs, this->introspection().component_indices.velocities[d]);
                    std::vector<bool> velocity_should_be_constrained(n_dofs);
                    std::vector<double> velocity_solution(n_dofs);
                    it->constrain_solution(cell,
                                           input_data.evaluation_points,
                                           velocity_component_indices,
                                           velocity_should_be_constrained,
                                           velocity_solution);

                    computed_quantities[q][2] = std::max(computed_quantities[q][2], double(velocity_should_be_constrained[q]));

                    if (velocity_should_be_constrained[q] == true)
                      computed_quantities[q][3+d] = velocity_solution[q]; // for values
                  }

                // for composition
                for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                  {
                    std::vector<unsigned int> composition_component_indices(n_dofs, this->introspection().component_indices.compositional_fields[c]);
                    std::vector<bool> composition_should_be_constrained(n_dofs);
                    std::vector<double> composition_solution(n_dofs);
                    it->constrain_solution(cell,
                                           input_data.evaluation_points,
                                           composition_component_indices,
                                           composition_should_be_constrained,
                                           composition_solution);

                    computed_quantities[q][3+dim] = std::max(computed_quantities[q][3+dim], double(composition_should_be_constrained[q]));

                    if (composition_should_be_constrained[q] == true)
                      computed_quantities[q][4+dim+c] = composition_solution[q]; // for values
                  }
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
                                                  "A visualization output object that visualizes "
                                                  "the should-be-constrained indicators "
                                                  "of the prescribed solution.")

    }
  }
}
