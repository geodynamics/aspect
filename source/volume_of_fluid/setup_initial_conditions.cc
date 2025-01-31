/*
 Copyright (C) 2016 - 2020 by the authors of the ASPECT code.

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

#include <aspect/global.h>
#include <aspect/simulator.h>
#include <aspect/volume_of_fluid/handler.h>
#include <aspect/volume_of_fluid/utilities.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/base/parsed_function.h>

namespace aspect
{
  template <int dim>
  void VolumeOfFluidHandler<dim>::set_initial_volume_fractions ()
  {
    for (unsigned int f=0; f<n_volume_of_fluid_fields; ++f)
      {

        switch (initialization_data_type[f])
          {
            case VolumeOfFluid::VolumeOfFluidInputType::composition:
              initialize_from_composition_field (data[f]);
              break;
            case VolumeOfFluid::VolumeOfFluidInputType::level_set:
              initialize_from_level_set (data[f]);
              break;
            default:
              Assert(false, ExcNotImplemented ());
          }

        const unsigned int volume_of_fluidN_blockidx = data[f].reconstruction.block_index;
        const unsigned int volume_of_fluidLS_blockidx = data[f].level_set.block_index;
        update_volume_of_fluid_normals (data[f], sim.solution);
        sim.old_solution.block(volume_of_fluidN_blockidx) = sim.solution.block(volume_of_fluidN_blockidx);
        sim.old_old_solution.block(volume_of_fluidN_blockidx) = sim.solution.block(volume_of_fluidN_blockidx);
        sim.old_solution.block(volume_of_fluidLS_blockidx) = sim.solution.block(volume_of_fluidLS_blockidx);
        sim.old_old_solution.block(volume_of_fluidLS_blockidx) = sim.solution.block(volume_of_fluidLS_blockidx);

        // Update associated composition field
        const typename Simulator<dim>::AdvectionField composition_field = Simulator<dim>::AdvectionField::composition(data[f].composition_index);
        update_volume_of_fluid_composition (composition_field, data[f], sim.solution);
        const unsigned int volume_of_fluid_C_blockidx = composition_field.block_index(this->introspection());
        sim.old_solution.block(volume_of_fluid_C_blockidx) = sim.solution.block(volume_of_fluid_C_blockidx);
        sim.old_old_solution.block(volume_of_fluid_C_blockidx) = sim.solution.block(volume_of_fluid_C_blockidx);
      }
  }

  template <int dim>
  void VolumeOfFluidHandler<dim>::initialize_from_composition_field (
    const VolumeOfFluidField<dim> &field)
  {
    LinearAlgebra::BlockVector initial_solution;

    initial_solution.reinit(sim.system_rhs, false);

    const QIterated<dim> quadrature (QMidpoint<1>(), n_init_samples);
    FEValues<dim, dim> fe_init (this->get_mapping(), this->get_fe(), quadrature,
                                update_JxW_values | update_quadrature_points);

    std::vector<types::global_dof_index>
    local_dof_indices (this->get_fe().dofs_per_cell);

    const FEVariable<dim> &volume_of_fluid_var = field.volume_fraction;
    const unsigned int component_index = volume_of_fluid_var.first_component_index;
    const unsigned int blockidx = volume_of_fluid_var.block_index;
    const unsigned int volume_of_fluid_ind
      = this->get_fe().component_to_system_index(component_index, 0);

    // Initialize state based on provided function
    for (auto cell : this->get_dof_handler().active_cell_iterators ())
      {
        if (!cell->is_locally_owned ())
          continue;

        // Calculate approximation for volume
        cell->get_dof_indices (local_dof_indices);

        fe_init.reinit (cell);

        double volume_of_fluid_val = 0.0;
        double cell_vol = 0.0;

        for (unsigned int q = 0; q < fe_init.n_quadrature_points; ++q)
          {
            const double fraction_at_point = this->get_initial_composition_manager().initial_composition(fe_init.quadrature_point(q),
                                             field.composition_index);
            const double JxW = fe_init.JxW(q);
            volume_of_fluid_val += fraction_at_point * JxW;
            cell_vol += JxW;
          }

        volume_of_fluid_val /= cell_vol;

        volume_of_fluid_val = std::min(volume_of_fluid_val, 1.0);
        volume_of_fluid_val = std::max(volume_of_fluid_val, 0.0);

        initial_solution (local_dof_indices[volume_of_fluid_ind]) = volume_of_fluid_val;
      }

    initial_solution.compress(VectorOperation::insert);

    // Apply constraints and update solution blocks. This is duplicated from
    // Simulator<dim>::set_initial_temperature_and_compositional_fields()
    // in order to separate the volume of fluid algorithm as much as possible from
    // the rest of the code.
    sim.compute_current_constraints();
    sim.current_constraints.distribute(initial_solution);

    sim.solution.block(blockidx) = initial_solution.block(blockidx);
    sim.old_solution.block(blockidx) = initial_solution.block(blockidx);
    sim.old_old_solution.block(blockidx) = initial_solution.block(blockidx);
  }

  template <int dim>
  void VolumeOfFluidHandler<dim>::initialize_from_level_set (
    const VolumeOfFluidField<dim> &field)
  {
    LinearAlgebra::BlockVector initial_solution;

    initial_solution.reinit(sim.system_rhs, false);

    const QIterated<dim> quadrature (QMidpoint<1>(), n_init_samples);
    FEValues<dim, dim> fe_init (this->get_mapping(),
                                this->get_fe(),
                                quadrature,
                                update_JxW_values | update_quadrature_points);

    const double h = 1.0/n_init_samples;

    std::vector<types::global_dof_index>
    local_dof_indices (this->get_fe().dofs_per_cell);

    const FEVariable<dim> &volume_of_fluid_var = field.volume_fraction;
    const unsigned int component_index = volume_of_fluid_var.first_component_index;
    const unsigned int blockidx = volume_of_fluid_var.block_index;
    const unsigned int volume_of_fluid_ind
      = this->get_fe().component_to_system_index(component_index, 0);

    // Initialize state based on provided function
    for (auto cell : this->get_dof_handler().active_cell_iterators ())
      {
        if (!cell->is_locally_owned ())
          continue;

        // Calculate approximation for volume
        cell->get_dof_indices (local_dof_indices);

        const double cell_diam = cell->diameter();
        const double d_func = this->get_initial_composition_manager().initial_composition(cell->barycenter(),
                              field.composition_index);
        fe_init.reinit (cell);

        double volume_of_fluid_val = 0.0;
        double cell_vol = 0.0;

        if (d_func <=-0.5*cell_diam)
          {
            volume_of_fluid_val = 0.0;
          }
        else if (d_func >= 0.5*cell_diam)
          {
            volume_of_fluid_val = 1.0;
          }
        else
          {

            // For each quadrature point compute an approximation to the fluid fraction in the surrounding region
            for (unsigned int q = 0; q < fe_init.n_quadrature_points; ++q)
              {
                double d = 0.0;
                Tensor<1, dim, double> grad;
                Point<dim> xU = quadrature.point (q);

                // Get an approximation to local normal at the closest interface (level set gradient)
                // and the distance to the closest interface (value of level set function)
                for (unsigned int di = 0; di < dim; ++di)
                  {
                    Point<dim> xH, xL;
                    xH = xU;
                    xL = xU;
                    xH[di] += 0.5*h;
                    xL[di] -= 0.5*h;
                    const double dH = this->get_initial_composition_manager().initial_composition(cell->intermediate_point(xH),
                                                                                                  field.composition_index);
                    const double dL = this->get_initial_composition_manager().initial_composition(cell->intermediate_point(xL),
                                                                                                  field.composition_index);
                    grad[di] = (dL-dH);
                    d += (0.5/dim)*(dH+dL);
                  }
                // Use the basic fluid fraction formula to compute an approximation to the fluid fraction
                const double fraction_at_point = VolumeOfFluid::Utilities::compute_fluid_fraction (grad, d);
                const double JxW = fe_init.JxW(q);
                volume_of_fluid_val += fraction_at_point * JxW;
                cell_vol += JxW;
              }
            volume_of_fluid_val /= cell_vol;
          }

        initial_solution (local_dof_indices[volume_of_fluid_ind]) = volume_of_fluid_val;
      }

    initial_solution.compress(VectorOperation::insert);

    // Apply constraints and update solution blocks. This is duplicated from
    // Simulator<dim>::set_initial_temperature_and_compositional_fields()
    // in order to separate the volume of fluid algorithm as much as possible from
    // the rest of the code.
    sim.compute_current_constraints();
    sim.current_constraints.distribute(initial_solution);

    sim.solution.block(blockidx) = initial_solution.block(blockidx);
    sim.old_solution.block(blockidx) = initial_solution.block(blockidx);
    sim.old_old_solution.block(blockidx) = initial_solution.block(blockidx);
  }
}

namespace aspect
{
#define INSTANTIATE(dim) \
  template void VolumeOfFluidHandler<dim>::set_initial_volume_fractions ();\
  template void VolumeOfFluidHandler<dim>::initialize_from_composition_field (const VolumeOfFluidField<dim> &field); \
  template void VolumeOfFluidHandler<dim>::initialize_from_level_set (const VolumeOfFluidField<dim> &field);

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
