/*
  Copyright (C) 2011, 2012, 2013, 2014 by the authors of the ASPECT code.

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
/*  $Id$  */

#include <aspect/simulator.h>
#include <aspect/global.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>


using namespace dealii;


namespace aspect
{

  template <int dim>
  void Simulator<dim>::free_surface_execute()
  {
    if (!parameters.free_surface_enabled)
      return;
    pcout << "FS: free_surface_execute()" << std::endl;








    free_surface_displace_mesh();
  }

  template <int dim>
  void Simulator<dim>::free_surface_make_constraints()
  {
    if (!parameters.free_surface_enabled)
      return;
    pcout << "FS: free_surface_make_constraints()" << std::endl;

    mesh_constraints.clear();

//Assert(we need a free surface)

  }

  template <int dim>
  void Simulator<dim>::free_surface_setup_dofs()
  {
    if (!parameters.free_surface_enabled)
      return;
    pcout << "FS: free_surface_setup_dofs()" << std::endl;

    // these live in the same FE as the velocity variable:
    mesh_velocity.reinit(introspection.index_sets.system_partitioning, introspection.index_sets.system_relevant_partitioning, mpi_communicator);
    old_mesh_velocity.reinit(introspection.index_sets.system_partitioning, introspection.index_sets.system_relevant_partitioning, mpi_communicator);


    free_surface_dof_handler.distribute_dofs(free_surface_fe);

    pcout << "FS: n_dofs = " << free_surface_dof_handler.n_dofs() << std::endl;

    // Renumber the DoFs hierarchical so that we get the
    // same numbering if we resume the computation. This
    // is because the numbering depends on the order the
    // cells are created.
    DoFRenumbering::hierarchical (free_surface_dof_handler);
//  DoFRenumbering::component_wise (free_surface_dof_handler,
//      introspection.components_to_blocks);



    mesh_locally_owned = free_surface_dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs (free_surface_dof_handler,
                                             mesh_locally_relevant);

    mesh_vertices.reinit(mesh_locally_owned, mesh_locally_relevant, mpi_communicator);
    //if we are just starting, we need to initialize mesh_vertices
    if (this->timestep_number == 0)
      {
        pcout << "FS: get initial mesh vertices" << std::endl;

        LinearAlgebra::Vector distributed_mesh_vertices;
        distributed_mesh_vertices.reinit(mesh_locally_owned, mpi_communicator);

        const std::vector<Point<dim> > mesh_support_points
          = free_surface_fe.base_element(0).get_unit_support_points();
        FEValues<dim> mesh_points (mapping, free_surface_fe,
                                   mesh_support_points, update_quadrature_points);
        std::vector<unsigned int> cell_dof_indices (free_surface_fe.dofs_per_cell);

        typename DoFHandler<dim>::active_cell_iterator cell = free_surface_dof_handler.begin_active(),
                                                       endc = free_surface_dof_handler.end();
        for (; cell != endc; ++cell)
          if (cell->is_locally_owned())
            {
              mesh_points.reinit(cell);
              cell->get_dof_indices (cell_dof_indices);
              for (unsigned int j=0; j<free_surface_fe.base_element(0).dofs_per_cell; ++j)
                for (unsigned int dir=0; dir<dim; ++dir)
                  {
                    unsigned int support_point_index
                      = free_surface_fe.component_to_system_index(/*velocity component=*/ dir,
                                                                                          /*dof index within component=*/ j);
                    distributed_mesh_vertices[cell_dof_indices[support_point_index]] = mesh_points.quadrature_point(j)[dir];
                  }
            }

        distributed_mesh_vertices.compress(VectorOperation::insert);
        mesh_vertices = distributed_mesh_vertices;
      }

    free_surface_make_constraints();

    // matrix
    {
      mesh_matrix.clear ();

      Table<2,DoFTools::Coupling> coupling (dim, dim);
      coupling.fill(DoFTools::none);

      for (unsigned int c=0; c<dim; ++c)
        coupling[c][c] = DoFTools::always;

#ifdef USE_PETSC
      LinearAlgebra::CompressedSparsityPattern sp(error);

#else
      TrilinosWrappers::SparsityPattern sp (mesh_locally_owned,mesh_locally_owned,
                                            mpi_communicator);
#endif

      DoFTools::make_sparsity_pattern (free_surface_dof_handler,
                                       coupling, sp,
                                       mesh_constraints, false,
                                       Utilities::MPI::
                                       this_mpi_process(mpi_communicator));

#ifdef USE_PETSC
      SparsityTools::distribute_sparsity_pattern(sp,
                                                 free_surface_dof_handler.locally_owned_dofs_per_processor(),
                                                 mpi_communicator, mesh_locally_relevant);

      sp.compress();

      mesh_matrix.reinit (mesh_locally_owned, mesh_locally_owned, sp, mpi_communicator);
#else
      sp.compress();

      mesh_matrix.reinit (sp);
#endif

    }


  }

  template <int dim>
  void Simulator<dim>::free_surface_displace_mesh()
  {
    pcout << "FS: free_surface_displace_mesh()" << std::endl;

    typename DoFHandler<dim>::active_cell_iterator  cell = free_surface_dof_handler.begin_active(),
                                                    endc = free_surface_dof_handler.end();

    for (cell = free_surface_dof_handler.begin_active(); cell != endc; ++cell)
      if (cell->is_artificial() == false)
        for (unsigned int vertex_no = 0; vertex_no < GeometryInfo<dim>::vertices_per_cell; ++vertex_no)
          {
            Point<dim> &v=cell->vertex(vertex_no);
            for (unsigned int dir=0; dir<dim; ++dir)
              v(dir) = mesh_vertices(
                         cell->vertex_dof_index(vertex_no, dir)
                       ); //enforce the vertex position
          }

  }

}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template class Simulator<dim>;

  ASPECT_INSTANTIATE(INSTANTIATE)
}
