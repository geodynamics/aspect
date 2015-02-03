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


#include <aspect/postprocess/pressure_statistics.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    PressureStatistics<dim>::execute (TableHandler &statistics)
    {
      // create a quadrature formula based on the pressure element alone.
      // we need to compute max and min of the pressure as well, which
      // may be on the boundary of the cell, so we use an iterated
      // trapezoidal rule instead of the usual Gauss rule
      //
      // iterate it 'degree' times to make sure our evaluation points are
      // in fact the support points.
      const QIterated<dim> quadrature_formula (QTrapez<1>(),
                                               this->get_fe().base_element(this->introspection().base_elements.pressure).degree);
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);

      std::vector<double> pressure_values(n_q_points);

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      double local_pressure_integral = 0;
      double local_min_pressure      = std::numeric_limits<double>::max();
      double local_max_pressure      = -std::numeric_limits<double>::max();

      // compute the integral quantities by quadrature. note that compared to
      // the temperature statistics postprocessor, we can not just loop over
      // the pressure DoFs because they may be intermingled with the
      // velocity DoFs if we use a direct solver
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            fe_values[this->introspection().extractors.pressure].get_function_values (this->get_solution(),
                                                                                      pressure_values);
            for (unsigned int q=0; q<n_q_points; ++q)
              {
                const double value = pressure_values[q];

                local_pressure_integral += value*fe_values.JxW(q);
                local_min_pressure = std::min (local_min_pressure, value);
                local_max_pressure = std::max (local_max_pressure, value);
              }
          }

      const double global_pressure_integral
        = Utilities::MPI::sum (local_pressure_integral, this->get_mpi_communicator());

      // now do the reductions that are
      // min/max operations. do them in
      // one communication by multiplying
      // one value by -1
      double global_min_pressure = 0;
      double global_max_pressure = 0;
      {
        double local_values[2] = { -local_min_pressure, local_max_pressure };
        double global_values[2];

        Utilities::MPI::max (local_values, this->get_mpi_communicator(), global_values);

        global_min_pressure = -global_values[0];
        global_max_pressure = global_values[1];
      }

      double global_mean_pressure = global_pressure_integral / this->get_volume();
      statistics.add_value ("Minimal pressure (Pa)",
                            global_min_pressure);
      statistics.add_value ("Average pressure (Pa)",
                            global_mean_pressure);
      statistics.add_value ("Maximal pressure (Pa)",
                            global_max_pressure);

      // also make sure that the other columns filled by the this object
      // all show up with sufficient accuracy and in scientific notation
      {
        const char *columns[] = { "Minimal pressure (Pa)",
                                  "Average pressure (Pa)",
                                  "Maximal pressure (Pa)"
                                };
        for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
          {
            statistics.set_precision (columns[i], 8);
            statistics.set_scientific (columns[i], true);
          }
      }

      std::ostringstream output;
      output.precision(4);
      output << global_min_pressure << " Pa, "
             << global_pressure_integral / this->get_volume() << " Pa, "
             << global_max_pressure << " Pa";

      return std::pair<std::string, std::string> ("Pressure min/avg/max:",
                                                  output.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(PressureStatistics,
                                  "pressure statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the pressure field.")
  }
}
