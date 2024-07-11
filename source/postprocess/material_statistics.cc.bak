/*
  Copyright (C) 2018 - 2022 by the authors of the ASPECT code.

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


#include <aspect/postprocess/material_statistics.h>
#include <aspect/material_model/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    MaterialStatistics<dim>::execute (TableHandler &statistics)
    {
      // create a quadrature formula based on the temperature element alone.
      const Quadrature<dim> &quadrature_formula = this->introspection().quadratures.temperature;
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
      in.requested_properties = MaterialModel::MaterialProperties::density | MaterialModel::MaterialProperties::viscosity;

      double local_volume = 0.0;
      double local_mass = 0.0;
      double local_viscosity = 0.0;

      // compute the integral quantities by quadrature
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            in.reinit(fe_values, cell, this->introspection(), this->get_solution());

            this->get_material_model().fill_additional_material_model_inputs(in, this->get_solution(), fe_values, this->introspection());
            this->get_material_model().evaluate(in, out);

            for (unsigned int q=0; q<n_q_points; ++q)
              {
                local_mass += out.densities[q] * fe_values.JxW(q);
                local_viscosity += out.viscosities[q] * fe_values.JxW(q);
                local_volume += fe_values.JxW(q);
              }
          }

      // compute the sum over all processors
      const double global_mass = Utilities::MPI::sum (local_mass, this->get_mpi_communicator());
      const double global_viscosity = Utilities::MPI::sum (local_viscosity, this->get_mpi_communicator());
      const double global_volume = Utilities::MPI::sum (local_volume, this->get_mpi_communicator());
      const double average_density = global_mass / global_volume;
      const double average_viscosity = global_viscosity / global_volume;

      const std::string name1("Average density (kg/m^3)");
      statistics.add_value (name1, average_density);
      statistics.set_precision (name1, 8);
      statistics.set_scientific (name1, true);

      const std::string name2("Average viscosity (Pa s)");
      statistics.add_value (name2, average_viscosity);
      statistics.set_precision (name2, 8);
      statistics.set_scientific (name2, true);

      const std::string name3("Total mass (kg)");
      statistics.add_value (name3, global_mass);
      statistics.set_precision (name3, 8);
      statistics.set_scientific (name3, true);

      std::ostringstream output;
      output.precision(4);

      output << average_density << " kg/m^3, "
             << average_viscosity << " Pa s, "
             << global_mass << " kg";

      return std::pair<std::string, std::string> ("Average density / Average viscosity / Total mass: ",
                                                  output.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(MaterialStatistics,
                                  "material statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the material properties. In particular, it computes the "
                                  "volume-averages of the density and viscosity, and the "
                                  "total mass in the model. Specifically, it implements "
                                  "the following formulas in each time step: "
                                  "$\\left<\\rho\\right> = \\frac{1}{|\\Omega|} \\int_\\Omega \\rho(\\mathbf x) \\, \\text{d}x$, "
                                  "$\\left<\\eta\\right> = \\frac{1}{|\\Omega|} \\int_\\Omega \\eta(\\mathbf x) \\, \\text{d}x$, "
                                  "$M = \\int_\\Omega \\rho(\\mathbf x) \\, \\text{d}x$, "
                                  "where $|\\Omega|$ is the volume of the domain.")
  }
}
