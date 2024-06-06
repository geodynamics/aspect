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


#include "isotherm_depth.h"
#include <aspect/geometry_model/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    IsothermDepth<dim>::execute (TableHandler &statistics)
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

      std::vector<Point<dim>> position_values(n_q_points);
      std::vector<double> temperature_values(n_q_points);

      double local_isotherm_depth = 0.0;
      const double max_depth = this->get_geometry_model().maximal_depth();
      const double dT = 1.1;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {

            fe_values.reinit (cell);

            position_values = fe_values.get_quadrature_points();

            fe_values[this->introspection().extractors.temperature].get_function_values (this->get_solution(),
                                                                                         temperature_values);

            for (unsigned int i=0; i<n_q_points; ++i)
              {
                // Calculate the biggest depth per processor
                if (temperature_values[i] > (isotherm_value - dT) && temperature_values[i] < (isotherm_value + dT) &&
                    (max_depth - position_values[i][dim-1]) > local_isotherm_depth)
                  {
                    local_isotherm_depth = max_depth - position_values[i][dim-1];
                  }
              }
          }

      // compute the maximum depth over all processors
      const  double isotherm_depth =
        Utilities::MPI::max (local_isotherm_depth, this->get_mpi_communicator());

      statistics.add_value ("Isotherm depth [m]", isotherm_depth);

      // also make sure that the other columns filled by this object
      // all show up with sufficient accuracy and in scientific notation
      {
        const char *columns[] = { "Isotherm depth [m]"
                                };
        for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
          {
            statistics.set_precision (columns[i], 8);
            statistics.set_scientific (columns[i], true);
          }
      }

      std::ostringstream output;
      output.precision(4);
      output << isotherm_depth
             << " m ";
      return std::pair<std::string, std::string> ("Isotherm depth [m]:",
                                                  output.str());
    }



    template <int dim>
    void IsothermDepth<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Isotherm depth");
        {
          prm.declare_entry("Isotherm value", "1573",
                            Patterns::Double(0),
                            "The temperature value of the isotherm "
                            "which maximum depth is tracked. Unit: [K].");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void IsothermDepth<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Isotherm depth");
        {
          isotherm_value = prm.get_double("Isotherm value");
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
    ASPECT_REGISTER_POSTPROCESSOR(IsothermDepth,
                                  "isotherm depth",
                                  "A postprocessor that computes the maximum depth "
                                  "of a user-given isotherm.")
  }
}
