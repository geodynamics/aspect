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
      AssertThrow(this->introspection().compositional_name_exists("porosity"),
                  ExcMessage("The 'fluid velocity statistics' postprocessor requires a "
                             "compositional field called porosity."));

      const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

      AssertThrow(this->get_parameters().include_melt_transport ||
                  this->get_parameters().compositional_field_methods[porosity_idx] == Parameters<dim>::AdvectionFieldMethod::fem_darcy_field,
                  ExcMessage("The 'fluid velocity statistics' postprocessor requires advecting the 'porosity' compositional "
                             "field by either the darcy velocity or the melt velocity."));

      const bool output_melt_velocity = this->get_parameters().include_melt_transport;

      // In this postprocessor, we only calculate the max fluid velocity. For this, we use a trapezoidal quadrature
      // which contains points on the borders of the cell.
      const QIterated<dim> quadrature_formula (QTrapezoid<1>(),
                                               this->get_parameters().stokes_velocity_degree);
      const unsigned int n_q_points = quadrature_formula.size();

      const UpdateFlags update_flags
        = UpdateFlags(
            !output_melt_velocity
            ?
            update_values |
            update_gradients |
            update_quadrature_points |
            update_JxW_values
            :
            update_values |
            update_JxW_values);

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_flags);

      std::vector<Tensor<1,dim>> solid_velocity_values(n_q_points);
      std::vector<Tensor<1,dim>> fluid_velocity_values(n_q_points);

      MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points, this->n_compositional_fields());

      MeltHandler<dim>::create_material_model_outputs(out);

      double local_max_fluid_velocity = 0;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);

            // If the porosity is advected with the melt velocity, use the melt velocity...
            if (output_melt_velocity)
              {
                const FEValuesExtractors::Vector ex_u_f = this->introspection().variable("fluid velocity").extractor_vector();
                fe_values[ex_u_f].get_function_values (this->get_solution(), fluid_velocity_values);
                for (unsigned int q = 0; q < n_q_points; ++q)
                  local_max_fluid_velocity = std::max (fluid_velocity_values[q].norm(),
                                                       local_max_fluid_velocity);
              }

            // Otherwise, use the Darcy velocity
            else
              {
                fe_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                                                                                            solid_velocity_values);
                in.reinit(fe_values, cell, this->introspection(), this->get_solution());

                this->get_material_model().evaluate(in, out);

                const std::shared_ptr<MaterialModel::MeltOutputs<dim>> fluid_out
                  = out.template get_additional_output_object<MaterialModel::MeltOutputs<dim>>();

                for (unsigned int q = 0; q < n_q_points; ++q)
                  {
                    const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(fe_values.quadrature_point(q));
                    const double porosity = std::max(in.composition[q][porosity_idx], 1e-10);
                    const double solid_density = out.densities[q];
                    const double fluid_density = fluid_out->fluid_densities[q];
                    const double fluid_viscosity = fluid_out->fluid_viscosities[q];
                    const double permeability = fluid_out->permeabilities[q];
                    const Tensor<1,dim> darcy_velocity = solid_velocity_values[q] -
                                                         (permeability / fluid_viscosity / porosity) *
                                                         gravity * (solid_density - fluid_density);
                    local_max_fluid_velocity = std::max (darcy_velocity.norm(),
                                                         local_max_fluid_velocity);
                  }
              }
          }

      const double global_max_fluid_velocity
        = Utilities::MPI::max (local_max_fluid_velocity, this->get_mpi_communicator());

      const std::string units = (this->convert_output_to_years() == true) ? "m/year" : "m/s";
      const double unit_scale_factor = (this->convert_output_to_years() == true) ? year_in_seconds : 1.0;
      const std::string advection_method = output_melt_velocity ? "melt" : "Darcy";
      const std::vector<std::string> column_name = {"Max. " + advection_method + " velocity (" + units + ")"
                                                   };

      statistics.add_value (column_name[0],
                            global_max_fluid_velocity * unit_scale_factor);

      // also make sure that the other columns filled by this object
      // all show up with sufficient accuracy and in scientific notation
      for (auto &column : column_name)
        {
          statistics.set_precision (column, 8);
          statistics.set_scientific (column, true);
        }

      std::ostringstream output;
      output.precision(3);
      output << global_max_fluid_velocity *unit_scale_factor
             << ' ' << units;

      return std::pair<std::string, std::string> ("Max " + advection_method + " velocity:",
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
                                  "A postprocessor that computes the maximum fluid velocity "
                                  "in the computational domain. The fluid velocity is either "
                                  "the melt velocity or the Darcy velocity depending on the "
                                  "advection method chosen for the porosity compositional field.")
  }
}
