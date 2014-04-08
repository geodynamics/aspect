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
/*  $Id: main.cc 2430 2014-04-08 14:41:18Z heister $  */

#include <aspect/simulator.h>
#include <aspect/global.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>


using namespace dealii;


namespace aspect
{

template <int dim>
void Simulator<dim>::free_surface_execute()
{
  if (!parameters.free_surface_enabled)
    return;
  pcout << "FS: free_surface_execute()" << std::endl;








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

}
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template class Simulator<dim>;

  ASPECT_INSTANTIATE(INSTANTIATE)
}
