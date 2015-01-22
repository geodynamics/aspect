/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include <aspect/postprocess/velocity_statistics.h>
#include <aspect/material_model/simple.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    VelocityStatistics<dim>::execute (TableHandler &statistics)
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

      double local_velocity_square_integral = 0;
      double local_max_velocity = 0;

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            fe_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                                                                                        velocity_values);
            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                local_velocity_square_integral += ((velocity_values[q] * velocity_values[q]) *
                                                   fe_values.JxW(q));
                local_max_velocity = std::max (std::sqrt(velocity_values[q]*velocity_values[q]),
                                               local_max_velocity);
              }
          }

      const double global_velocity_square_integral
        = Utilities::MPI::sum (local_velocity_square_integral, this->get_mpi_communicator());
      const double global_max_velocity
        = Utilities::MPI::max (local_max_velocity, this->get_mpi_communicator());

      const double vrms = std::sqrt(global_velocity_square_integral) /
                          std::sqrt(this->get_volume());

      if (this->convert_output_to_years() == true)
        {
          statistics.add_value ("RMS velocity (m/year)",
                                vrms * year_in_seconds);
          statistics.add_value ("Max. velocity (m/year)",
                                global_max_velocity * year_in_seconds);

          // also make sure that the other columns filled by the this object
          // all show up with sufficient accuracy and in scientific notation
          {
            const char *columns[] = { "RMS velocity (m/year)",
                                      "Max. velocity (m/year)"
                                    };
            for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
              {
                statistics.set_precision (columns[i], 8);
                statistics.set_scientific (columns[i], true);
              }
          }
        }
      else
        {
          statistics.add_value ("RMS velocity (m/s)", vrms);
          statistics.add_value ("Max. velocity (m/s)", global_max_velocity);

          // also make sure that the other columns filled by the this object
          // all show up with sufficient accuracy and in scientific notation
          {
            const char *columns[] = { "RMS velocity (m/s)",
                                      "Max. velocity (m/s)"
                                    };
            for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
              {
                statistics.set_precision (columns[i], 8);
                statistics.set_scientific (columns[i], true);
              }
          }
        }

      std::ostringstream output;
      output.precision(3);
      if (this->convert_output_to_years() == true)
        output << vrms *year_in_seconds
               << " m/year, "
               << global_max_velocity *year_in_seconds
               << " m/year";
      else
        output << vrms
               << " m/s, "
               << global_max_velocity
               << " m/s";

      return std::pair<std::string, std::string> ("RMS, max velocity:",
                                                  output.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(VelocityStatistics,
                                  "velocity statistics",
                                  "A postprocessor that computes some statistics about the "
                                  "velocity field.")
  }
}
