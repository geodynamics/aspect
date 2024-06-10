/*
  Copyright (C) 2024 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include "composition_trms_statistics.h"
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
    CompositionTrmsStatistics<dim>::execute (TableHandler &statistics)
    {
      // create a quadrature formula based on the temperature element alone.
      const Quadrature<dim> &quadrature_formula = this->introspection().quadratures.temperature;
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);

      std::vector<double> temperature_values(n_q_points);
      std::vector<double> compositional_values(n_q_points);

      std::vector<double> local_temperature_square_integral (this->n_compositional_fields());
      std::vector<double> local_area_integral (this->n_compositional_fields());
      double local_total_temperature_square_integral = 0.;
      double local_total_area_integral = 0.;
      const double kelvin_to_celsius = 273.15;
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);

            fe_values[this->introspection().extractors.temperature].get_function_values(this->get_solution(),
                                                                                        temperature_values);

            // Loop over the quadrature points to get the temperature and area
            // integral for the whole domain.
            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                local_total_temperature_square_integral += (((temperature_values[q]-kelvin_to_celsius) *
                                                             (temperature_values[q]-kelvin_to_celsius)) * fe_values.JxW(q));
                local_total_area_integral += fe_values.JxW(q);
              }

            // For each field, loop over the quadrature points to get the temperature
            // and area per field.
            for (unsigned int c = 0; c < this->n_compositional_fields(); ++c)
              {
                fe_values[this->introspection().extractors.compositional_fields[c]].get_function_values(this->get_solution(),
                    compositional_values);

                for (unsigned int q = 0; q < n_q_points; ++q)
                  {
                    if (compositional_values[q] >= 0.5)
                      {
                        local_temperature_square_integral[c] += (((temperature_values[q]-kelvin_to_celsius) *
                                                                  (temperature_values[q]-kelvin_to_celsius)) * fe_values.JxW(q));
                        local_area_integral[c] += fe_values.JxW(q);
                      }
                  }
              }
          }

      std::vector<double> global_temperature_square_integral (local_temperature_square_integral.size());
      Utilities::MPI::sum (local_temperature_square_integral, this->get_mpi_communicator(), global_temperature_square_integral);
      std::vector<double> global_area_integral (local_area_integral.size());
      Utilities::MPI::sum (local_area_integral, this->get_mpi_communicator(), global_area_integral);
      const double global_total_temperature_square_integral = Utilities::MPI::sum(local_total_temperature_square_integral, this->get_mpi_communicator());
      const double global_total_area_integral = Utilities::MPI::sum(local_total_area_integral, this->get_mpi_communicator());

      // compute the RMS temperature for each compositional field and for the selected compositional fields combined
      std::vector<double> Trms_per_composition(local_area_integral.size());
      double temperature_square_integral_selected_fields = 0., area_integral_selected_fields = 0.;

      for (unsigned int c = 0; c < this->n_compositional_fields(); ++c)
        {
          Trms_per_composition[c] = std::sqrt(global_temperature_square_integral[c]) /
                                    std::sqrt(global_area_integral[c]);

          const std::vector<std::string>::iterator selected_field_it = std::find(selected_fields.begin(), selected_fields.end(), this->introspection().name_for_compositional_index(c));
          if (selected_field_it != selected_fields.end())
            {
              temperature_square_integral_selected_fields += global_temperature_square_integral[c];
              area_integral_selected_fields += global_area_integral[c];
            }
        }

      const double Trms_selected_fields = std::sqrt(temperature_square_integral_selected_fields) / std::sqrt(area_integral_selected_fields);

      // compute the RMS temperature over the whole domain
      const double Trms_whole_domain = std::sqrt(global_total_temperature_square_integral) / std::sqrt(global_total_area_integral);

      // finally produce something for the statistics file
      for (unsigned int c = 0; c < this->n_compositional_fields(); ++c)
        {
          statistics.add_value("RMS temperature (C) for composition " + this->introspection().name_for_compositional_index(c),
                               Trms_per_composition[c]);

          // also make sure that the other columns filled by this object
          // all show up with sufficient accuracy and in scientific notation
          const std::string columns[] = {"RMS temperature (C) for composition " + this->introspection().name_for_compositional_index(c)};

          for (unsigned int i = 0; i < sizeof(columns) / sizeof(columns[0]); ++i)
            {
              statistics.set_precision(columns[i], 8);
              statistics.set_scientific(columns[i], true);
            }
        }

      // Also output the RMS temperature for the selected fields
      statistics.add_value("RMS temperature (C) for the selected fields",
                           Trms_selected_fields);

      const std::string column = {"RMS temperature (C) for the selected fields"};

      statistics.set_precision(column, 8);
      statistics.set_scientific(column, true);

      // Also output the RMS temperature for the whole domain
      statistics.add_value("RMS temperature (C) for the whole domain",
                           Trms_whole_domain);

      const std::string column_whole_domain = {"RMS temperature (C) for the whole domain"};

      statistics.set_precision(column_whole_domain, 8);
      statistics.set_scientific(column_whole_domain, true);

      std::ostringstream output;
      output.precision(4);

      for (unsigned int c = 0; c < this->n_compositional_fields(); ++c)
        {
          output << Trms_per_composition[c]
                 << " C";
          output << " // ";
        }
      output << Trms_selected_fields << " C // ";
      output << Trms_whole_domain << " C";

      return std::pair<std::string, std::string>("RMS temperature for compositions, combined selected fields and whole domain:",
                                                 output.str());
    }



    template <int dim>
    void CompositionTrmsStatistics<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Composition RMS temperature statistics");
        {
          prm.declare_entry("Names of selected compositional fields", "",
                            Patterns::List(Patterns::Anything()),
                            "A list of names for each of the compositional fields that "
                            "you want to compute the combined RMS temperature for.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void CompositionTrmsStatistics<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Composition RMS temperature statistics");
        {
          selected_fields = Utilities::split_string_list(prm.get("Names of selected compositional fields"));

          AssertThrow((selected_fields.size() > 0) &&
                      (selected_fields.size() <= this->n_compositional_fields()),
                      ExcMessage("The length of the list of names for the compositional "
                                 "fields for which the RMS temperature is to be summed must be larger than zero "
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
    ASPECT_REGISTER_POSTPROCESSOR(CompositionTrmsStatistics,
                                  "composition RMS temperature statistics",
                                  "A postprocessor that computes the root-mean-square (RMS) temperature "
                                  "over the area spanned by each compositional field (i.e. where "
                                  "the field values are larger or equal to 0.5). It also computes the RMS "
                                  "temperature over the summed area of a list of fields given by the user, "
                                  "and over the whole domain.")
  }
}
