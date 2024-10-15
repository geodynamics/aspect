/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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


#include <aspect/mesh_refinement/volume_of_fluid_interface.h>
#include <aspect/volume_of_fluid/handler.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/geometry_info.h>

#include <deal.II/numerics/error_estimator.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_tools.h>

namespace aspect
{
  namespace MeshRefinement
  {

    template <int dim>
    void
    VolumeOfFluidInterface<dim>::tag_additional_cells() const
    {
      // Break early if DoFs have not been distributed yet.
      if (this->get_dof_handler().n_dofs() == 0)
        return;

      const QMidpoint<dim> qMidC;

      // Create a map from vertices to adjacent cells
      const std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>
      vertex_to_cells(GridTools::vertex_to_cell_map(this->get_triangulation()));

      std::set<typename Triangulation<dim>::active_cell_iterator> marked_cells;
      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               qMidC,
                               update_values |
                               update_quadrature_points);

      // Create block vector to indicate when other mpi process believes cell
      // borders an interface cell and therefore neighbors should be refined
      // Use the first vof block, due to being the correct size, and only
      // needing one indicator
      LinearAlgebra::BlockVector interface_contained_local(
        this->introspection().index_sets.system_partitioning,
        this->get_mpi_communicator());

      const FiniteElement<dim> &system_fe = this->get_fe();

      const VolumeOfFluidField<dim> vof_field = this->get_volume_of_fluid_handler().field_struct_for_field_index(0);
      const unsigned int volume_fraction_block = vof_field.volume_fraction.block_index;
      const unsigned int volume_of_fluid_c_index = vof_field.volume_fraction.first_component_index;
      const unsigned int volume_of_fluid_ind
        = this->get_fe().component_to_system_index(volume_of_fluid_c_index, 0);

      interface_contained_local.block(volume_fraction_block) = 0.0;

      std::vector<types::global_dof_index> local_dof_indices (system_fe.dofs_per_cell);

      for (unsigned int f=0; f<this->get_volume_of_fluid_handler().get_n_fields(); ++f)
        {

          const double volume_fraction_threshold = this->get_volume_of_fluid_handler().get_volume_fraction_threshold();

          const FEValuesExtractors::Scalar volume_of_fluid_field = this->get_volume_of_fluid_handler().field_struct_for_field_index(f)
                                                                   .volume_fraction.extractor_scalar();

          std::vector<double> volume_of_fluid_values(qMidC.size());
          std::vector<double> neighbor_volume_of_fluid_values(qMidC.size());


          for (const auto &cell : this->get_dof_handler().active_cell_iterators())
            if (!cell->is_artificial())
              {
                bool refine_current_cell = false;

                // Get cell volume_of_fluid
                fe_values.reinit(cell);
                fe_values[volume_of_fluid_field].get_function_values(this->get_solution(),
                                                                     volume_of_fluid_values);

                // Handle overshoots
                volume_of_fluid_values[0] = std::min(volume_of_fluid_values[0], 1.0);
                volume_of_fluid_values[0] = std::max(volume_of_fluid_values[0], 0.0);

                // Check if at interface
                if (volume_of_fluid_values[0] > volume_fraction_threshold && volume_of_fluid_values[0] < (1.0 - volume_fraction_threshold))
                  {
                    refine_current_cell = true;
                  }

                if (!refine_current_cell)
                  {
                    for (const unsigned int f : cell->face_indices())
                      {
                        const bool cell_has_periodic_neighbor = cell->has_periodic_neighbor(f);
                        const typename DoFHandler<dim>::face_iterator face = cell->face(f);

                        // Skip if face is at boundary, and does not have a periodic neighbor
                        if (face->at_boundary() && !cell_has_periodic_neighbor)
                          continue;

                        const typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor_or_periodic_neighbor(f);

                        Assert(cell.state()==IteratorState::valid, ExcInternalError());

                        if ((!face->at_boundary() && !face->has_children())
                            ||
                            (face->at_boundary() && cell->periodic_neighbor_is_coarser(f))
                            ||
                            (face->at_boundary() && neighbor->level() == cell->level() && neighbor->is_active()))
                          {
                            if (neighbor->is_active() && !neighbor->is_artificial())
                              {
                                fe_values.reinit(neighbor);
                                fe_values[volume_of_fluid_field].get_function_values(this->get_solution(),
                                                                                     neighbor_volume_of_fluid_values);

                                // Handle overshoots
                                neighbor_volume_of_fluid_values[0] = std::min(neighbor_volume_of_fluid_values[0], 1.0);
                                neighbor_volume_of_fluid_values[0] = std::max(neighbor_volume_of_fluid_values[0], 0.0);

                                if (std::abs(neighbor_volume_of_fluid_values[0]-volume_of_fluid_values[0])>volume_fraction_threshold)
                                  {
                                    refine_current_cell = true;
                                    break;
                                  }
                              }
                          }
                        else
                          {
                            for (unsigned int subface=0; subface < face->n_children(); ++subface)
                              {
                                const typename DoFHandler<dim>::active_cell_iterator neighbor_sub =
                                  (cell_has_periodic_neighbor
                                   ?
                                   cell->periodic_neighbor_child_on_subface(f, subface)
                                   :
                                   cell->neighbor_child_on_subface(f, subface));

                                Assert(neighbor_sub.state()==IteratorState::valid, ExcInternalError());

                                if (neighbor_sub->is_artificial())
                                  continue;

                                fe_values.reinit(neighbor_sub);
                                fe_values[volume_of_fluid_field].get_function_values(this->get_solution(),
                                                                                     neighbor_volume_of_fluid_values);

                                // Handle overshoots
                                neighbor_volume_of_fluid_values[0] = std::min(neighbor_volume_of_fluid_values[0], 1.0);
                                neighbor_volume_of_fluid_values[0] = std::max(neighbor_volume_of_fluid_values[0], 0.0);

                                if (std::abs(neighbor_volume_of_fluid_values[0]-volume_of_fluid_values[0])>volume_fraction_threshold)
                                  {
                                    refine_current_cell = true;
                                    break;
                                  }

                              }
                            if (refine_current_cell)
                              break;
                          }
                      }
                  }

                if (refine_current_cell)
                  {
                    // Fractional volume
                    marked_cells.insert(cell);
                    if (cell->is_locally_owned())
                      {
                        cell->clear_coarsen_flag ();
                        cell->set_refine_flag ();
                        // Mark in vector, will be true here if true on any process
                        cell->get_dof_indices(local_dof_indices);
                        interface_contained_local(local_dof_indices[volume_of_fluid_ind]) = 1.0;
                      }
                  }

              }
        }

      // Now communicate and mark any cells not already included, this could be
      // reduced to only loop over cells bordering another process
      LinearAlgebra::BlockVector interface_contained_global(
        this->introspection().index_sets.system_partitioning,
        this->introspection().index_sets.system_relevant_partitioning,
        this->get_mpi_communicator());

      interface_contained_global.block(volume_fraction_block) = interface_contained_local.block(volume_fraction_block);
      interface_contained_global.update_ghost_values();

      const FEValuesExtractors::Scalar ic_extract = this->get_volume_of_fluid_handler().field_struct_for_field_index(0)
                                                    .volume_fraction.extractor_scalar();
      std::vector<double> ic_values(qMidC.size());

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (!cell->is_artificial())
          {
            fe_values.reinit(cell);
            fe_values[ic_extract].get_function_values(interface_contained_global,
                                                      ic_values);

            // If the cell has been indicated to need refinement add it
            if (ic_values[0]>0.5)
              {
                marked_cells.insert(cell);
              }
          }

      // Now mark for refinement all cells that are a neighbor of a cell that contains the interface

      std::set<typename Triangulation<dim>::active_cell_iterator> marked_cells_and_neighbors = marked_cells;
      typename std::set<typename parallel::distributed::Triangulation<dim>::active_cell_iterator>::const_iterator mcells = marked_cells.begin(),
                                                                                                                  endmc = marked_cells.end();
      for (; mcells != endmc; ++mcells)
        {
          typename parallel::distributed::Triangulation<dim>::active_cell_iterator mcell = *mcells;
          for (const unsigned int vertex_index : mcell->vertex_indices())
            {
              const std::set<typename Triangulation<dim>::active_cell_iterator> &neighbor_cells = vertex_to_cells[mcell->vertex_index(vertex_index)];
              for (const auto &neighbor_cell : neighbor_cells)
                {
                  if (neighbor_cell->is_active() && neighbor_cell->is_locally_owned())
                    {
                      neighbor_cell->clear_coarsen_flag ();
                      neighbor_cell->set_refine_flag ();
                      marked_cells_and_neighbors.insert(neighbor_cell);
                    }
                }
            }

          // Check for periodic neighbors, and refine if existing
          for (const unsigned int f : mcell->face_indices())
            {
              if (mcell->has_periodic_neighbor(f))
                {
                  typename Triangulation<dim>::cell_iterator itr_tmp = mcell->periodic_neighbor(f);

                  if (itr_tmp->is_active() && itr_tmp->is_locally_owned())
                    {
                      itr_tmp->clear_coarsen_flag ();
                      itr_tmp->set_refine_flag ();
                      marked_cells_and_neighbors.insert(itr_tmp);
                    }
                }
            }
        }

      if (strict_coarsening)
        {
          for (const auto &cell : this->get_dof_handler().active_cell_iterators())
            if (cell->is_locally_owned())
              {
                if (marked_cells_and_neighbors.find(cell) != marked_cells_and_neighbors.end())
                  {
                    //Refinement already requested
                  }
                else
                  {
                    if (cell->is_active())
                      {
                        cell->set_coarsen_flag();
                        cell->clear_refine_flag();
                      }
                  }
              }

        }

    }

    template <int dim>
    void
    VolumeOfFluidInterface<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Volume of fluid interface");
        {
          prm.declare_entry("Strict coarsening", "false",
                            Patterns::Bool(),
                            "If true, then explicitly coarsen any cells not "
                            "neighboring the VolumeOfFluid interface.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    VolumeOfFluidInterface<dim>::parse_parameters (ParameterHandler &prm)
    {
      AssertThrow(this->get_parameters().volume_of_fluid_tracking_enabled,
                  ExcMessage("The 'volume_of_fluid boundary' mesh refinement strategy requires that the 'Use interface tracking' parameter is enabled."));

      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Volume of fluid interface");
        {
          strict_coarsening = prm.get_bool("Strict coarsening");
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
  namespace MeshRefinement
  {
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(VolumeOfFluidInterface,
                                              "volume of fluid interface",
                                              "A class that implements a mesh refinement criterion, which "
                                              "ensures a minimum level of refinement near the volume of fluid interface boundary.")
  }
}
