/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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


#include <aspect/postprocess/composition_velocity_statistics.h>
#include <aspect/global.h>
#include <aspect/utilities.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    CompositionVelocityStatistics<dim>::execute (TableHandler &statistics)
    {
      const QGauss<dim> quadrature_formula (this->get_fe()
                                            .base_element(this->introspection().base_elements.velocities).degree+1);
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);

      std::vector<Tensor<1,dim>> velocity_values(n_q_points);
      std::vector<double> compositional_values(n_q_points);

      std::vector<double> local_velocity_square_integral(this->n_compositional_fields());
      std::vector<double> local_area_integral(this->n_compositional_fields());

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            fe_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                                                                                        velocity_values);

            for (unsigned int c = 0; c < this->n_compositional_fields(); ++c)
              {
                fe_values[this->introspection().extractors.compositional_fields[c]].get_function_values(this->get_solution(),
                    compositional_values);

                for (unsigned int q = 0; q < n_q_points; ++q)
                  {
                    if (compositional_values[q] >= 0.5)
                      {
                        local_velocity_square_integral[c] += ((velocity_values[q] * velocity_values[q]) *
                                                              fe_values.JxW(q));
                        local_area_integral[c] += fe_values.JxW(q);
                      }
                  }
              }
          }

      std::vector<double> global_velocity_square_integral(local_velocity_square_integral.size());
      Utilities::MPI::sum(local_velocity_square_integral, this->get_mpi_communicator(), global_velocity_square_integral);
      std::vector<double> global_area_integral(local_area_integral.size());
      Utilities::MPI::sum(local_area_integral, this->get_mpi_communicator(), global_area_integral);

      // compute the RMS velocity for each compositional field and for the selected compositonal fields combined
      std::vector<double> vrms_per_composition(local_area_integral.size());
      double velocity_square_integral_selected_fields = 0., area_integral_selected_fields = 0.;
      for (unsigned int c = 0; c < this->n_compositional_fields(); ++c)
        {
          vrms_per_composition[c] = std::sqrt(global_velocity_square_integral[c]) /
                                    std::sqrt(global_area_integral[c]);

          const std::vector<std::string>::iterator selected_field_it = std::find(selected_fields.begin(), selected_fields.end(), this->introspection().name_for_compositional_index(c));
          if (selected_field_it != selected_fields.end())
            {
              velocity_square_integral_selected_fields += global_velocity_square_integral[c];
              area_integral_selected_fields += global_area_integral[c];
            }
        }

      const double vrms_selected_fields = std::sqrt(velocity_square_integral_selected_fields) / std::sqrt(area_integral_selected_fields);

      const std::string unit = (this->convert_output_to_years()) ? "m/year" : "m/s";
      const double time_scaling = (this->convert_output_to_years()) ? year_in_seconds : 1.0;

      // finally produce something for the statistics file
      for (unsigned int c = 0; c < this->n_compositional_fields(); ++c)
        {
          statistics.add_value("RMS velocity (" + unit + ") for composition " + this->introspection().name_for_compositional_index(c),
                               time_scaling * vrms_per_composition[c]);

          // also make sure that the other columns filled by this object
          // all show up with sufficient accuracy and in scientific notation
          const std::string columns[] = {"RMS velocity (" + unit + ") for composition " +
                                         this->introspection().name_for_compositional_index(c)
                                        };

          for (const auto &column : columns)
            {
              statistics.set_precision(column, 8);
              statistics.set_scientific(column, true);
            }
        }

      // Also output the selected fields vrms
      statistics.add_value("RMS velocity (" + unit + ") for the selected field ",
                           time_scaling * vrms_selected_fields);

      const std::string column = {"RMS velocity (" + unit + ") for the selected field "};

      statistics.set_precision(column, 8);
      statistics.set_scientific(column, true);

      std::ostringstream output;
      output.precision(4);

      for (unsigned int c = 0; c < this->n_compositional_fields(); ++c)
        {
          output << time_scaling *vrms_per_composition[c]
                 << " " << unit;
          output << " // ";
        }
      output << time_scaling *vrms_selected_fields;

      return std::pair<std::string, std::string>("RMS velocity for compositions and combined selected fields:",
                                                 output.str());
    }



    template <int dim>
    void
    CompositionVelocityStatistics<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Composition velocity statistics");
        {
          prm.declare_entry("Names of selected compositional fields", "",
                            Patterns::List(Patterns::Anything()),
                            "A list of names for each of the compositional fields that "
                            "you want to compute the combined RMS velocity for.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    CompositionVelocityStatistics<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Composition velocity statistics");
        {
          selected_fields = Utilities::split_string_list(prm.get("Names of selected compositional fields"));

          AssertThrow((selected_fields.size() > 0) &&
                      (selected_fields.size() <= this->n_compositional_fields()),
                      ExcMessage("The length of the list of names for the compositional "
                                 "fields for which the RMS velocity is to be summed must be larger than zero "
                                 "and smaller or equal to the number of compositional fields."));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(CompositionVelocityStatistics,
                                  "composition velocity statistics",
                                  "A postprocessor that computes the root mean square velocity "
                                  "over the area spanned by each compositional field (i.e. where "
                                  "the field values are larger or equal to 0.5.")
  }
}
