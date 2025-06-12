/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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


#include <aspect/postprocess/fluid_velocity_statistics.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/global.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    FluidVelocityStatistics<dim>::execute (TableHandler &statistics)
    {
      const Quadrature<dim> &quadrature_formula = this->introspection().quadratures.velocities;
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_gradients |
                               update_quadrature_points |
                               update_JxW_values);
      std::vector<Tensor<1,dim>> velocity_values(n_q_points);
      std::vector<Tensor<1,dim>> fluid_velocity_values(n_q_points);

      MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points, this->n_compositional_fields());
      MeltHandler<dim>::create_material_model_outputs(out);
      MaterialModel::MeltOutputs<dim> *fluid_out = out.template get_additional_output<MaterialModel::MeltOutputs<dim>>();

      double local_fluid_velocity_square_integral = 0;
      double local_max_fluid_velocity = 0;

      unsigned int porosity_idx = numbers::invalid_unsigned_int;
      if (this->introspection().composition_type_exists(CompositionalFieldDescription::porosity))
        {
          porosity_idx = this->introspection().find_composition_type(CompositionalFieldDescription::porosity);
        }

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            fe_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                                                                                        velocity_values);
            in.reinit(fe_values, cell, this->introspection(), this->get_solution());
            this->get_material_model().evaluate(in, out);
            fluid_out = out.template get_additional_output<MaterialModel::MeltOutputs<dim>>();
            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(fe_values.quadrature_point(q));
                const double porosity = std::max(in.composition[q][porosity_idx], 1e-10);
                const double solid_density = out.densities[q];
                const double fluid_density = fluid_out->fluid_densities[q];
                const double fluid_viscosity = fluid_out->fluid_viscosities[q];
                const double permeability = fluid_out->permeabilities[q];
                fluid_velocity_values[q] = velocity_values[q] -
                                           (permeability / fluid_viscosity / porosity) *
                                           gravity * (solid_density - fluid_density);
                local_fluid_velocity_square_integral += ((fluid_velocity_values[q] * fluid_velocity_values[q]) *
                                                         fe_values.JxW(q));
                local_max_fluid_velocity = std::max (std::sqrt(fluid_velocity_values[q]*fluid_velocity_values[q]),
                                                     local_max_fluid_velocity);
              }
          }

      const double global_fluid_velocity_square_integral
        = Utilities::MPI::sum (local_fluid_velocity_square_integral, this->get_mpi_communicator());
      const double global_max_fluid_velocity
        = Utilities::MPI::max (local_max_fluid_velocity, this->get_mpi_communicator());

      const double fluid_v_rms = std::sqrt(global_fluid_velocity_square_integral) /
                                 std::sqrt(this->get_volume());

      const std::string units = (this->convert_output_to_years() == true) ? "m/year" : "m/s";
      const double unit_scale_factor = (this->convert_output_to_years() == true) ? year_in_seconds : 1.0;
      const std::vector<std::string> column_names = {"RMS fluid velocity (" + units + ")",
                                                     "Max. fluid velocity (" + units + ")"
                                                    };

      statistics.add_value (column_names[0],
                            fluid_v_rms * unit_scale_factor);
      statistics.add_value (column_names[1],
                            global_max_fluid_velocity * unit_scale_factor);

      // also make sure that the other columns filled by this object
      // all show up with sufficient accuracy and in scientific notation
      for (auto &column : column_names)
        {
          statistics.set_precision (column, 8);
          statistics.set_scientific (column, true);
        }

      std::ostringstream output;
      output.precision(3);
      output << fluid_v_rms *unit_scale_factor
             << ' ' << units << ", "
             << global_max_fluid_velocity *unit_scale_factor
             << ' ' << units;

      return std::pair<std::string, std::string> ("RMS fluid, max fluid velocity:",
                                                  output.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(FluidVelocityStatistics,
                                  "fluid velocity statistics",
                                  "A postprocessor that computes the root mean square and "
                                  "maximum velocity in the computational domain.")
  }
}
