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



#include <aspect/postprocess/melt_statistics.h>
#include <aspect/melt.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    MeltStatistics<dim>::execute (TableHandler &statistics)
    {

      const QIterated<dim> quadrature_formula_for_max (QTrapezoid<1>(),
                                                       this->get_parameters().temperature_degree);
      const unsigned int n_q_points_for_max = quadrature_formula_for_max.size();

      const Quadrature<dim> &quadrature_formula_for_integral = this->introspection().quadratures.temperature;
      const unsigned int n_q_points_for_integral = quadrature_formula_for_integral.size();

      FEValues<dim> fe_values_for_max (this->get_mapping(),
                                       this->get_fe(),
                                       quadrature_formula_for_max,
                                       update_values   |
                                       update_quadrature_points |
                                       update_JxW_values);

      FEValues<dim> fe_values_for_integral (this->get_mapping(),
                                            this->get_fe(),
                                            quadrature_formula_for_integral,
                                            update_values   |
                                            update_quadrature_points |
                                            update_JxW_values);

      MaterialModel::MaterialModelInputs<dim> in_for_max(fe_values_for_max.n_quadrature_points, this->n_compositional_fields());
      MaterialModel::MaterialModelInputs<dim> in_for_integral(fe_values_for_integral.n_quadrature_points, this->n_compositional_fields());


      std::ostringstream output;
      output.precision(4);

      double local_melt_integral = 0.0;
      double local_min_melt = std::numeric_limits<double>::max();
      double local_max_melt = std::numeric_limits<double>::lowest();

      // compute the integral quantities by quadrature
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            // fill material model inputs
            fe_values_for_max.reinit (cell);
            fe_values_for_integral.reinit (cell);

            in_for_max.reinit(fe_values_for_max, cell, this->introspection(), this->get_solution());
            in_for_integral.reinit(fe_values_for_max, cell, this->introspection(), this->get_solution());

            // we can only postprocess melt fractions if the material model that is used
            // in the simulation has implemented them
            // otherwise, set them to zero
            std::vector<double> melt_fractions_for_max(n_q_points_for_max, 0.0);
            std::vector<double> melt_fractions_for_integral(n_q_points_for_integral, 0.0);

            if (MaterialModel::MeltFractionModel<dim>::is_melt_fraction_model(this->get_material_model()))
              {
                MaterialModel::MeltFractionModel<dim>::as_melt_fraction_model(this->get_material_model())
                .melt_fractions(in_for_max, melt_fractions_for_max);
                MaterialModel::MeltFractionModel<dim>::as_melt_fraction_model(this->get_material_model())
                .melt_fractions(in_for_integral, melt_fractions_for_integral);
              }

            for (unsigned int q=0; q<n_q_points_for_max; ++q)
              {
                local_min_melt       = std::min(local_min_melt, melt_fractions_for_max[q]);
                local_max_melt       = std::max(local_max_melt, melt_fractions_for_max[q]);
              }
            for (unsigned int q=0; q<n_q_points_for_integral; ++q)
              {
                local_melt_integral += melt_fractions_for_integral[q] * fe_values_for_integral.JxW(q);
              }

          }

      const double global_melt_integral
        = Utilities::MPI::sum (local_melt_integral, this->get_mpi_communicator());
      double global_min_melt = 0;
      double global_max_melt = 0;

      // now do the reductions that are
      // min/max operations. do them in
      // one communication by multiplying
      // one value by -1
      {
        double local_values[2] = { -local_min_melt, local_max_melt };
        double global_values[2];

        Utilities::MPI::max (local_values, this->get_mpi_communicator(), global_values);

        global_min_melt = -global_values[0];
        global_max_melt = global_values[1];
      }


      // finally produce something for the statistics file
      statistics.add_value ("Minimal melt fraction",
                            global_min_melt);
      statistics.add_value ("Total melt volume",
                            global_melt_integral);
      statistics.add_value ("Maximal melt fraction",
                            global_max_melt);

      // also make sure that the other columns filled by this object
      // all show up with sufficient accuracy and in scientific notation
      {
        const char *columns[] = { "Minimal melt fraction",
                                  "Total melt volume",
                                  "Maximal melt fraction"
                                };
        for (auto &column : columns)
          {
            statistics.set_precision (column, 8);
            statistics.set_scientific (column, true);
          }
      }

      output << global_min_melt << ", "
             << global_melt_integral << ", "
             << global_max_melt;

      return std::pair<std::string, std::string> ("Melt fraction min/total/max:",
                                                  output.str());

    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(MeltStatistics,
                                  "melt statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the melt (volume) fraction. If the material model "
                                  "does not implement a melt fraction function, the output is "
                                  "set to zero.")
  }
}
