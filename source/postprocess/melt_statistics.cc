/*
  Copyright (C) 2015 - 2016 by the authors of the ASPECT code.

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

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    MeltStatistics<dim>::execute (TableHandler &statistics)
    {
      // create a quadrature formula based on the temperature element alone.
      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.temperature).degree+1);
      const unsigned int n_q_points = quadrature_formula.size();
      std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),
                                                            std::vector<double> (n_q_points));

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_gradients |
                               update_quadrature_points |
                               update_JxW_values);

      MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, this->n_compositional_fields());

      std::ostringstream output;
      output.precision(4);

      double local_melt_integral = 0.0;
      double local_min_melt = std::numeric_limits<double>::max();
      double local_max_melt = -std::numeric_limits<double>::max();

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      // compute the integral quantities by quadrature
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            // fill material model inputs
            fe_values.reinit (cell);
            in.reinit(fe_values, cell, this->introspection(), this->get_solution());

            // we can only postprocess melt fractions if the material model that is used
            // in the simulation has implemented them
            // otherwise, set them to zero
            std::vector<double> melt_fractions(n_q_points, 0.0);
            if (const MaterialModel::MeltFractionModel<dim> *
                melt_material_model = dynamic_cast <const MaterialModel::MeltFractionModel<dim>*> (&this->get_material_model()))
              melt_material_model->melt_fractions(in, melt_fractions);

            for (unsigned int q=0; q<n_q_points; ++q)
              {
                local_melt_integral += melt_fractions[q] * fe_values.JxW(q);
                local_min_melt       = std::min(local_min_melt, melt_fractions[q]);
                local_max_melt       = std::max(local_max_melt, melt_fractions[q]);
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
        for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
          {
            statistics.set_precision (columns[i], 8);
            statistics.set_scientific (columns[i], true);
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
