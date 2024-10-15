/*
  Copyright (C) 2021 - 2024 by the authors of the ASPECT code.

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


#include <aspect/postprocess/max_depth_field.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    MaxDepthField<dim>::execute (TableHandler &statistics)
    {
      // create a quadrature formula based on the compositional element alone.
      // be defensive about determining that a compositional field actually exists
      AssertThrow(this->n_compositional_fields() > 0,
                  ExcMessage("This postprocessor cannot be used without compositional fields."));
      const Quadrature<dim> &quadrature_formula = this->introspection().quadratures.compositional_field_max;

      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points);

      std::vector<double> compositional_values(n_q_points);
      std::vector<Point<dim>> position_values(n_q_points);

      std::vector<double> local_max_depth(this->n_compositional_fields(), std::numeric_limits<double>::lowest());

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);

            position_values = fe_values.get_quadrature_points();

            for (unsigned int c = 0; c < this->n_compositional_fields(); ++c)
              {
                fe_values[this->introspection().extractors.compositional_fields[c]].get_function_values(this->get_solution(), compositional_values);

                for (unsigned int q = 0; q < n_q_points; ++q)
                  {
                    const double depth = this->get_geometry_model().depth(position_values[q]);
                    if (compositional_values[q] >= 0.5 && depth > local_max_depth[c])
                      local_max_depth[c] = depth;
                  }
              }
          }

      // compute the max over all processors
      std::vector<double> global_max_depth(this->n_compositional_fields());
      Utilities::MPI::max (local_max_depth,
                           this->get_mpi_communicator(),
                           global_max_depth);

      // finally produce something for the statistics file
      for (unsigned int c = 0; c < this->n_compositional_fields(); ++c)
        {
          statistics.add_value("Max depth [m] for composition " + this->introspection().name_for_compositional_index(c),
                               global_max_depth[c]);
        }

      // also make sure that the other columns filled by this object
      // all show up with sufficient accuracy and in scientific notation
      for (unsigned int c = 0; c < this->n_compositional_fields(); ++c)
        {
          const std::string columns[] = {"Max depth [m] for composition " + this->introspection().name_for_compositional_index(c)};
          for (const auto &col: columns)
            {
              statistics.set_precision(col, 8);
              statistics.set_scientific(col, true);
            }
        }

      std::ostringstream output;
      output.precision(4);
      for (unsigned int c = 0; c < this->n_compositional_fields(); ++c)
        {
          output << global_max_depth[c];
          if (c + 1 != this->n_compositional_fields())
            output << " // ";
        }
      return std::pair<std::string, std::string>("Compositions max depth [m]:",
                                                 output.str());

    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(MaxDepthField,
                                  "maximum depth of field",
                                  "A postprocessor that for each compositional field "
                                  "outputs the largest depth at which "
                                  "a quadrature point is found where the field "
                                  "has a value of 0.5 or larger. For fields that do not "
                                  "represent materials, but for example track a certain quantity "
                                  "like strain, this value of 0.5 does not necessarily make sense. ")
  }
}
