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


#include "composition_velocity_statistics.h"
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

      std::vector<Tensor<1,dim> > velocity_values(n_q_points);
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

      // compute the RMS velocity for each compositional field and for the slab
      std::vector<double> vrms_per_composition(local_area_integral.size());
      double velocity_square_integral_slab = 0., area_integral_slab = 0.;
      for (unsigned int c = 0; c < this->n_compositional_fields(); ++c)
        {
          vrms_per_composition[c] = std::sqrt(global_velocity_square_integral[c]) /
                                    std::sqrt(global_area_integral[c]);

          const std::vector<std::string>::iterator slab_it = std::find(slab_compositions.begin(), slab_compositions.end(), this->introspection().name_for_compositional_index(c));
          if (slab_it != slab_compositions.end())
            {
              std::cout << "slab compo " << *slab_it << std::endl;
              velocity_square_integral_slab += global_velocity_square_integral[c];
              area_integral_slab += global_area_integral[c];
            }
        }

      const double vrms_slab = std::sqrt(velocity_square_integral_slab) / std::sqrt(area_integral_slab);

      std::string unit = (this->convert_output_to_years()) ? "m/year" : "m/s";

      // finally produce something for the statistics file
      for (unsigned int c = 0; c < this->n_compositional_fields(); ++c)
        {
          if (this->convert_output_to_years())
            statistics.add_value("RMS velocity (" + unit + ") for composition " + this->introspection().name_for_compositional_index(c),
                                 year_in_seconds * vrms_per_composition[c]);
          else
            statistics.add_value("RMS velocity (" + unit + ") for composition " + this->introspection().name_for_compositional_index(c),
                                 vrms_per_composition[c]);

          // also make sure that the other columns filled by this object
          // all show up with sufficient accuracy and in scientific notation
          const std::string columns[] = {"RMS velocity (" + unit + ") for composition " + this->introspection().name_for_compositional_index(c)};

          for (unsigned int i = 0; i < sizeof(columns) / sizeof(columns[0]); ++i)
            {
              statistics.set_precision(columns[i], 8);
              statistics.set_scientific(columns[i], true);
            }
        }

      // Also output the slab vrms
      if (this->convert_output_to_years())
        statistics.add_value("RMS velocity (" + unit + ") for slab ",
                             year_in_seconds * vrms_slab);
      else
        statistics.add_value("RMS velocity (" + unit + ") for slab ",
                             vrms_slab);

      const std::string column = {"RMS velocity (" + unit + ") for slab "};

      statistics.set_precision(column, 8);
      statistics.set_scientific(column, true);

      std::ostringstream output;
      output.precision(4);

      for (unsigned int c = 0; c < this->n_compositional_fields(); ++c)
        {
          if (this->convert_output_to_years())
            output << year_in_seconds *vrms_per_composition[c]
                   << " " << unit;
          else
            output << vrms_per_composition[c]
                   << " " << unit;

          output << " // ";
        }
      if (this->convert_output_to_years())
        output << year_in_seconds *vrms_slab;
      else
        output << vrms_slab;

      return std::pair<std::string, std::string>("RMS velocity for compositions and slab:",
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
          prm.declare_entry("Names of slab compositional fields", "",
                            Patterns::List(Patterns::Anything()),
                            "A list of names for each of the compositional fields that. "
                            "makes up the subduction plate.");
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
          slab_compositions = Utilities::split_string_list(prm.get("Names of slab compositional fields"));

          AssertThrow((slab_compositions.size() > 0) &&
                      (slab_compositions.size() <= this->n_compositional_fields()),
                      ExcMessage("The length of the list of names for the compositional "
                                 "fields that make up the slabs much be larger than zero "
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
