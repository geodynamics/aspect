// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2015-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------



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
      // Use a Gauss-Lobatto quadrature formula based on the temperature
      // degree for computing the min/max, both of which may lie on the
      // boundaries of the cell.
      const QGaussLobatto<dim> quadrature_formula(this->get_fe().base_element(this->introspection().base_elements.temperature).degree + 2);
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_gradients |
                               update_quadrature_points |
                               update_JxW_values);

      MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points, this->n_compositional_fields());
      MeltHandler<dim>::create_material_model_outputs(out);

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
            fe_values.reinit (cell);

            in.reinit(fe_values, cell, this->introspection(), this->get_solution());
            this->get_material_model().evaluate(in, out);

            // we can only postprocess melt fractions if the material model that is used
            // in the simulation has implemented them
            // otherwise, set them to zero
            std::vector<double> melt_fractions(n_q_points, numbers::signaling_nan<double>());

            if (MaterialModel::MeltFractionModel<dim>::is_melt_fraction_model(this->get_material_model()))
              MaterialModel::MeltFractionModel<dim>::as_melt_fraction_model(this->get_material_model())
              .melt_fractions(in, melt_fractions, &out);

            for (unsigned int q=0; q<n_q_points; ++q)
              {
                local_min_melt       = std::min(local_min_melt, melt_fractions[q]);
                local_max_melt       = std::max(local_max_melt, melt_fractions[q]);
                local_melt_integral += melt_fractions[q] * fe_values.JxW(q);
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
