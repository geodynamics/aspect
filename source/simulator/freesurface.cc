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
#include <deal.II/numerics/vector_tools.h>


using namespace dealii;


namespace aspect
{

  template <int dim>
  void Simulator<dim>::free_surface_execute()
  {
    if (!parameters.free_surface_enabled)
      return;
    pcout << "FS: free_surface_execute()" << std::endl;

    free_surface_make_constraints();

    free_surface_solve_elliptic_problem();

    free_surface_calculate_mesh_displacement();

    free_surface_displace_mesh();
  }

  template <int dim>
  void Simulator<dim>::free_surface_make_constraints()
  {
    if (!parameters.free_surface_enabled)
      return;
    pcout << "FS: free_surface_make_constraints()" << std::endl;

    mesh_constraints.clear();
    mesh_constraints.reinit(mesh_locally_relevant);
    DoFTools::make_hanging_node_constraints(free_surface_dof_handler, mesh_constraints);

    //Add the vanilla periodic boundary constraints
    typedef std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int> > periodic_boundary_pairs;
    periodic_boundary_pairs pbp = geometry_model->get_periodic_boundary_pairs();
    for(periodic_boundary_pairs::iterator p = pbp.begin(); p != pbp.end(); ++p)
      DoFTools::make_periodicity_constraints(free_surface_dof_handler, (*p).first.first, (*p).first.second, (*p).second, mesh_constraints);


    //Zero out the displacement for the zero-velocity boundary indicators
    for (std::set<types::boundary_id>::const_iterator p = parameters.zero_velocity_boundary_indicators.begin();
         p != parameters.zero_velocity_boundary_indicators.end(); ++p)
      VectorTools::interpolate_boundary_values (free_surface_dof_handler, *p,
      ZeroFunction<dim>(dim), mesh_constraints);

    //make the tangential boundary indicators no displacement normal to the boundary
    VectorTools::compute_no_normal_flux_constraints (free_surface_dof_handler,
                                                     /* first_vector_component= */
                                                     0,
                                                     parameters.tangential_velocity_boundary_indicators,
                                                     mesh_constraints, mapping);

    //make the periodic boundary indicators no displacement normal to the boundary
    std::set< types::boundary_id > periodic_boundaries;
    for(periodic_boundary_pairs::iterator p = pbp.begin(); p != pbp.end(); ++p)
    {
      periodic_boundaries.insert((*p).first.first);
      periodic_boundaries.insert((*p).first.second);
    }
    VectorTools::compute_no_normal_flux_constraints (free_surface_dof_handler,
                                                     /* first_vector_component= */
                                                     0,
                                                     periodic_boundaries,
                                                     mesh_constraints, mapping);

    //For the free surface indicators we constrain the displacement to be v.n
    LinearAlgebra::Vector boundary_normal_velocity;
    boundary_normal_velocity.reinit(mesh_locally_owned, mesh_locally_relevant, mpi_communicator);
    free_surface_project_normal_velocity_onto_boundary( boundary_normal_velocity );

    //now insert the relevant part of the solution into the mesh constraints
    IndexSet constrained_dofs;
    DoFTools::extract_boundary_dofs(free_surface_dof_handler, ComponentMask(dim, true),
                                    constrained_dofs, parameters.free_surface_boundary_indicators);
    for( unsigned int i = 0; i < constrained_dofs.n_elements();  ++i)
    {
      types::global_dof_index index = constrained_dofs.nth_index_in_set(i);
      if (mesh_constraints.can_store_line(index))
        if(mesh_constraints.is_constrained(index)==false)
        {
          mesh_constraints.add_line(index);
          mesh_constraints.set_inhomogeneity(index, boundary_normal_velocity[index]);
        }
    }

    mesh_constraints.close();
  }


  template <int dim>
  void Simulator<dim>::free_surface_project_normal_velocity_onto_boundary(LinearAlgebra::Vector &output)
  {
    // TODO: should we use the extrapolated solution?
    pcout << "FS: free_surface_project_normal_velocity_onto_boundary()" << std::endl;



//    std::pair<double, Point<dim> > corrections = free_surface_determine_mesh_corrections();
//     double patch = corrections.first;
//     Point<dim> centroid = corrections.second;

     //stuff for iterating over the mesh
     QGauss<dim-1> face_quadrature(2);
     UpdateFlags update_flags = UpdateFlags(update_values | update_normal_vectors | update_JxW_values);
     FEFaceValues<dim> fs_fe_face_values (mapping, free_surface_fe, face_quadrature, update_flags);
     FEFaceValues<dim> fe_face_values (mapping, finite_element, face_quadrature, update_flags);
     const unsigned int n_face_q_points = fe_face_values.n_quadrature_points,
                        dofs_per_cell = fs_fe_face_values.dofs_per_cell;

     //stuff for assembling system
     std::vector<unsigned int> cell_dof_indices (dofs_per_cell);
     Vector<double> cell_vector (dofs_per_cell);
     FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

     //stuff for getting the velocity values
     std::vector<Tensor<1,dim> > velocity_values(n_face_q_points);

     //set up constraints
     ConstraintMatrix mass_matrix_constraints(mesh_locally_relevant);
     DoFTools::make_hanging_node_constraints(free_surface_dof_handler, mass_matrix_constraints);

     typedef std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int> > periodic_boundary_pairs;
     periodic_boundary_pairs pbp = geometry_model->get_periodic_boundary_pairs();
     for(periodic_boundary_pairs::iterator p = pbp.begin(); p != pbp.end(); ++p)
       DoFTools::make_periodicity_constraints(free_surface_dof_handler,
           (*p).first.first, (*p).first.second, (*p).second, mass_matrix_constraints);

     mass_matrix_constraints.close();

     //set up the matrix
     LinearAlgebra::SparseMatrix mass_matrix;
     TrilinosWrappers::SparsityPattern sp (mesh_locally_owned, mesh_locally_owned, mpi_communicator);
     DoFTools::make_sparsity_pattern (free_surface_dof_handler, sp, mass_matrix_constraints, false,
                                      Utilities::MPI::this_mpi_process(mpi_communicator));
     sp.compress();
     mass_matrix.reinit (sp);

     FEValuesExtractors::Vector extract_vel(0);

     //make distributed vectors.
     LinearAlgebra::Vector rhs, dist_solution;
     rhs.reinit(mesh_locally_owned, mpi_communicator);
     dist_solution.reinit(mesh_locally_owned, mpi_communicator);

     typename DoFHandler<dim>::active_cell_iterator
     cell = dof_handler.begin_active(), endc= dof_handler.end();
     typename DoFHandler<dim>::active_cell_iterator
     fscell = free_surface_dof_handler.begin_active();

     for (; cell!=endc; ++cell, ++fscell)
       if (cell->at_boundary() && cell->is_locally_owned())
         for(unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
           if( cell->face(face_no)->at_boundary() &&
             ((parameters.free_surface_boundary_indicators.find(cell->face(face_no)->boundary_indicator())
                != parameters.free_surface_boundary_indicators.end())))
           {
             fscell->get_dof_indices (cell_dof_indices);
             fs_fe_face_values.reinit (fscell, face_no);
             fe_face_values.reinit (cell, face_no);
             fe_face_values[introspection.extractors.velocities].get_function_values(solution, velocity_values);

             cell_vector = 0;
             cell_matrix = 0;
             for (unsigned int point=0; point<n_face_q_points; ++point)
               for (unsigned int i=0; i<dofs_per_cell; ++i)
               {
                 for (unsigned int j=0; j<dofs_per_cell; ++j)
                 {
                       cell_matrix(i,j) += (fs_fe_face_values[extract_vel].value(j,point) *
                           fs_fe_face_values[extract_vel].value(i,point) ) *
                           fs_fe_face_values.JxW(point);
                 }

                 cell_vector(i) += (fs_fe_face_values[extract_vel].value(i,point) *
                     fs_fe_face_values.normal_vector(point) ) *
                                   (velocity_values[point]*fs_fe_face_values.normal_vector(point)) *
                                   fs_fe_face_values.JxW(point);
               }

             mass_matrix_constraints.distribute_local_to_global (cell_matrix, cell_vector,
                                                          cell_dof_indices, mass_matrix, rhs, false);
           }

     rhs.compress (VectorOperation::add);
     mass_matrix.compress(VectorOperation::add);

     TrilinosWrappers::PreconditionJacobi preconditioner_mass;
     preconditioner_mass.initialize(mass_matrix);

     SolverControl solver_control(5*rhs.size(), 1e-7*rhs.l2_norm());
     SolverCG<LinearAlgebra::Vector> cg(solver_control);
     cg.solve (mass_matrix, dist_solution, rhs, preconditioner_mass);
     pcout << "\t\tsolved, its = " << solver_control.last_step() << std::endl;

     mass_matrix_constraints.distribute (dist_solution);
     output = dist_solution;
  }


  template <int dim>
  void Simulator<dim>::free_surface_solve_elliptic_problem()
  {
    pcout << "FS: free_surface_solve_elliptic_problem()" << std::endl;
    QGauss<dim> quadrature(1+1);
    UpdateFlags update_flags = UpdateFlags(update_values | update_JxW_values | update_gradients);
    FEValues<dim> fe_values (mapping, free_surface_fe, quadrature, update_flags);

    const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                       dofs_per_face = finite_element.dofs_per_face,
                       n_q_points    = fe_values.n_quadrature_points;

    std::vector<unsigned int> cell_dof_indices (dofs_per_cell);
    std::vector<unsigned int> face_dof_indices (dofs_per_face);
    Vector<double> cell_vector (dofs_per_cell);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    mesh_matrix = 0;

    //carry out the solution
    FEValuesExtractors::Vector extract_vel(0);

    LinearAlgebra::Vector rhs, poisson_solution;
    rhs.reinit(mesh_locally_owned, mpi_communicator);
    poisson_solution.reinit(mesh_locally_owned, mpi_communicator);

    typename DoFHandler<dim>::active_cell_iterator cell = free_surface_dof_handler.begin_active(),
        endc= free_surface_dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
      {
        cell->get_dof_indices (cell_dof_indices);
        fe_values.reinit (cell);

        cell_vector = 0;
        cell_matrix = 0;
        for (unsigned int point=0; point<n_q_points; ++point)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += scalar_product( fe_values[extract_vel].gradient(j,point),
                                    fe_values[extract_vel].gradient(i,point) ) *
                                    fe_values.JxW(point);
          }

        mesh_constraints.distribute_local_to_global (cell_matrix, cell_vector,
                                                     cell_dof_indices, mesh_matrix, rhs, false);
      }

    rhs.compress (VectorOperation::add);
    mesh_matrix.compress (VectorOperation::add);

    //Make the AMG preconditioner
    std::vector<std::vector<bool> > constant_modes;
    DoFTools::extract_constant_modes (free_surface_dof_handler,
                                      ComponentMask(dim, true),
                                      constant_modes);
    // TODO: think about keeping object between time steps
    LinearAlgebra::PreconditionAMG preconditioner_stiffness;
    LinearAlgebra::PreconditionAMG::AdditionalData Amg_data;
    Amg_data.constant_modes = constant_modes;
    Amg_data.elliptic = true;
    Amg_data.higher_order_elements = false;
    Amg_data.smoother_sweeps = 2;
    Amg_data.aggregation_threshold = 0.02;
    preconditioner_stiffness.initialize(mesh_matrix);

    SolverControl solver_control(5*rhs.size(), parameters.linear_stokes_solver_tolerance*rhs.l2_norm());
    SolverCG<LinearAlgebra::Vector> cg(solver_control);

    cg.solve (mesh_matrix, poisson_solution, rhs, preconditioner_stiffness);
    pcout << "\t\tsolved, its = " << solver_control.last_step() << std::endl;

    mesh_constraints.distribute (poisson_solution);
    mesh_vertex_velocity = poisson_solution;
  }


  template <int dim>
  void Simulator<dim>::free_surface_calculate_mesh_displacement()
  {
    LinearAlgebra::Vector distributed_mesh_vertices(mesh_locally_owned, mpi_communicator);
    LinearAlgebra::Vector distributed_mesh_vertex_velocity(mesh_locally_owned, mpi_communicator);

    distributed_mesh_vertices = mesh_vertices;
    distributed_mesh_vertex_velocity = mesh_vertex_velocity;

    pcout << "FS: mesh velocity: " << distributed_mesh_vertex_velocity.l2_norm() << std::endl;

    //actually do the ALE thing
    distributed_mesh_vertices.sadd(1.0, time_step, distributed_mesh_vertex_velocity);

    mesh_vertices = distributed_mesh_vertices;

    // TODO: calculate mesh_velocity
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
    mesh_vertex_velocity.reinit(mesh_locally_owned, mesh_locally_relevant, mpi_communicator);

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
