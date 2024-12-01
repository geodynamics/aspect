/*
  Copyright (C) 2019 - 2024 by the authors of the ASPECT code.

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


#include <aspect/postprocess/entropy_viscosity_statistics.h>
#include <aspect/global.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    EntropyViscosityStatistics<dim>::execute (TableHandler &statistics)
    {
      Vector<float> entropy_viscosity(this->get_triangulation().n_active_cells());
      this->get_artificial_viscosity(entropy_viscosity);

      // The entropy viscosity is cell-wise constant, so a simple midpoint quadrature
      // will be sufficient to integrate the value over cells.
      const QMidpoint<dim> quadrature_formula;
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_JxW_values);

      double local_maximum_viscosity = 0.0;
      double local_integrated_viscosity = 0.0;
      double local_volume = 0.0;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);

            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                const double entropy_viscosity_cell = entropy_viscosity[cell->active_cell_index()];
                local_maximum_viscosity = std::max(local_maximum_viscosity,entropy_viscosity_cell);
                local_volume += fe_values.JxW(q);
                local_integrated_viscosity += entropy_viscosity_cell * fe_values.JxW(q);
              }
          }

      const double global_integrated_viscosity
        = Utilities::MPI::sum (local_integrated_viscosity, this->get_mpi_communicator());
      const double global_integrated_volume
        = Utilities::MPI::sum (local_volume, this->get_mpi_communicator());
      const double global_maximum_viscosity
        = Utilities::MPI::max (local_maximum_viscosity, this->get_mpi_communicator());

      const double average_viscosity = global_integrated_viscosity / global_integrated_volume;

      statistics.add_value ("Max entropy viscosity (W/(m*K))",
                            global_maximum_viscosity);
      statistics.add_value ("Average entropy viscosity (W/(m*K))",
                            average_viscosity);

      // also make sure that the columns filled by this plugin
      // show up with sufficient accuracy and in scientific notation
      const char *columns[] = { "Max entropy viscosity (W/(m*K))",
                                "Average entropy viscosity (W/(m*K))"
                              };
      for (auto &column : columns)
        {
          statistics.set_precision (column, 8);
          statistics.set_scientific (column, true);
        }

      std::ostringstream output;
      output.precision(3);
      output << global_maximum_viscosity
             << " W/(m*K), "
             << average_viscosity
             << " W/(m*K)";


      return std::pair<std::string, std::string> ("Max / average entropy viscosity:",
                                                  output.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(EntropyViscosityStatistics,
                                  "entropy viscosity statistics",
                                  "A postprocessor that computes the maximum and volume averaged"
                                  "entropy viscosity stabilization for the temperature field.")
  }
}
