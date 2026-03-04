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

      // todo_post
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

        const aspect::PrescribedSolution::Manager<dim> &prescribed_solution_manager
          = Plugins::get_plugin_as_type<const aspect::PrescribedSolution::Manager<dim>>(this->get_prescribed_solution());

        for (auto &model_name : prescribed_solution_manager.get_active_plugin_names())
          {
            // for indicators
            solution_names.push_back("prescribed " + model_name + " indicator");
            std::replace(solution_names.back().begin(),solution_names.back().end(),' ', '_');

            // for values
            if (model_name == "velocity function")
              {
                for (unsigned int d=0; d<dim; ++d)
                  {
                    solution_names.push_back("prescribed " + model_name);
                    std::replace(solution_names.back().begin(),solution_names.back().end(),' ', '_');
                  }
              }
            else
              {
                solution_names.push_back("prescribed " + model_name);
                std::replace(solution_names.back().begin(),solution_names.back().end(),' ', '_');

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
        for (auto &model_name : model_names)
          {
            // for indicators
            interpretation.push_back (DataComponentInterpretation::component_is_scalar);

            // for values
            if (model_name == "velocity function")
              {
                for (unsigned int d=0; d<dim; ++d)
                  interpretation.push_back (DataComponentInterpretation::component_is_part_of_vector);
              }
            else
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

        const std::vector<std::string> &plugin_names = prescribed_solution_manager.get_active_plugin_names();
        const auto &plugin_objects = prescribed_solution_manager.get_active_plugins();

        const unsigned int n_dofs = input_data.evaluation_points.size();


        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            unsigned output_index = 0;
            for (unsigned int i=0; i<model_names.size(); ++i, ++output_index, ++output_index)
              {
                auto &model_name = model_names[i];
                const std::vector<std::string>::const_iterator plugin = std::find(plugin_names.begin(),
                                                                                  plugin_names.end(),
                                                                                  model_name);

                AssertThrow(plugin != plugin_names.end(),
                            ExcMessage("The prescribed solution postprocessor was asked for outputing from plugin "
                                       "with the name <" + model_name + ">, but no such plugin is activated as prescribed solution."));

                unsigned j = std::distance(plugin_names.begin(),plugin);
                auto it = plugin_objects.begin();
                std::advance(it, j);

                typename DoFHandler<dim>::active_cell_iterator cell;

                if (model_name == "temperature function")
                  {
                    std::vector<unsigned int> component_indices(n_dofs, this->introspection().component_indices.temperature);
                    std::vector<bool> should_be_constrained(n_dofs);
                    std::vector<double> solution(n_dofs);
                    (*it)->constrain_solution(cell,
                                              input_data.evaluation_points,
                                              component_indices,
                                              should_be_constrained,
                                              solution);


                    computed_quantities[q][output_index] = should_be_constrained[q]; // for indicators

                    computed_quantities[q][output_index+1] = solution[q]; // for values
                  }

                if (model_name == "velocity function")
                  {
                    computed_quantities[q][output_index] = 0.0; // for indicator, initialize as zero

                    for (unsigned int d=0; d<dim; ++d)
                      {
                        // specify the component
                        std::vector<unsigned int> component_indices(n_dofs, this->introspection().component_indices.velocities[d]);
                        std::vector<bool> should_be_constrained(n_dofs);
                        std::vector<double> solution(n_dofs);
                        (*it)->constrain_solution(cell,
                                                  input_data.evaluation_points,
                                                  component_indices,
                                                  should_be_constrained,
                                                  solution);

                        computed_quantities[q][output_index] = std::max(computed_quantities[q][output_index], double(should_be_constrained[q]));

                        computed_quantities[q][output_index+1+d] = solution[q]; // for values
                      }

                    for (unsigned int d=0; d<dim; ++d)
                      ++output_index;
                  }

              }

          }
      }

      template <int dim>
      void
      PrescribedSolutionPostprocessor<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Prescribed solution");
            {
              const std::string pattern_of_names
                = "temperature function|velocity function";

              prm.declare_entry("List of model names",
                                "temperature function",
                                Patterns::MultipleSelection(pattern_of_names),
                                "A comma-separated list of prescribed solution models that "
                                "will be written as both indicators and values"
                                "The following prescribed solution models are available for output:\n\n"
                                +
                                pattern_of_names);
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }



      template <int dim>
      void
      PrescribedSolutionPostprocessor<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Prescribed solution");
            {
              model_names = Utilities::split_string_list(prm.get ("List of model names"));
              AssertThrow(Utilities::has_unique_entries(model_names),
                          ExcMessage("The list of strings for the parameter "
                                     "'Postprocess/Visualization/Prescribed solution/List of model names' contains entries more than once. "
                                     "This is not allowed. Please check your parameter file."));
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
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
                                                  "PrescribedSolutionPostprocessor",
                                                  "A visualization output object that visualizes "
                                                  "the should-be-constrained indicators "
                                                  "of the prescribed solution.")

    }
  }
}
