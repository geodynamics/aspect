/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include <aspect/simulator.h>
#include <aspect/global.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/numerics/vector_tools.h>


using namespace dealii;


namespace aspect
{

  template <int dim>
  Simulator<dim>::FreeSurfaceHandler::FreeSurfaceHandler( Simulator<dim> &simulator,
                                                          ParameterHandler &prm)
    : sim(simulator),  //reference to the simulator that owns the FreeSurfaceHandler
      free_surface_fe (FE_Q<dim>(1),dim), //Q1 elements which describe the mesh geometry
      free_surface_dof_handler (sim.triangulation)
  {
    parse_parameters(prm);
  }

  template <int dim>
  void Simulator<dim>::FreeSurfaceHandler::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection ("Free surface");
    {
      prm.declare_entry("Free surface stabilization theta", "0.5",
                        Patterns::Double(0,1),
                        "Theta parameter described in Kaus et. al. 2010. "
                        "An unstabilized free surface can overshoot its "
                        "equilibrium position quite easily and generate "
                        "unphysical results.  One solution is to use a "
                        "quasi-implicit correction term to the forces near the "
                        "free surface.  This parameter describes how much "
                        "the free surface is stabilized with this term, "
                        "where zero is no stabilization, and one is fully "
                        "implicit.");
    }
    prm.leave_subsection ();
  }

  template <int dim>
  void Simulator<dim>::FreeSurfaceHandler::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection ("Free surface");
    {
      free_surface_theta = prm.get_double("Free surface stabilization theta");
    }
    prm.leave_subsection ();
  }



  template <int dim>
  void Simulator<dim>::FreeSurfaceHandler::execute()
  {
    if (!sim.parameters.free_surface_enabled)
      return;
    sim.computing_timer.enter_section("FreeSurface");

    //Make the constraints for the elliptic problem.  On the free surface, we
    //constrain mesh velocity to be v.n, on free slip it is constrainted to
    //be tangential, and on no slip boundaries it is zero.
    make_constraints();

    //Assemble and solve the vector Laplace problem which determines
    //the mesh velocitiy in the interior of the domain
    solve_elliptic_problem();

    //calculate the displacement from the mesh velocity based on the
    //current time step size.  Also interpolates the Q1 mesh velocity
    //system onto the Stokes velocity system for ALE corrections
    calculate_mesh_displacement();

    //Actually loop over the mesh vertices and displace them.
    //Also called after redistribution with p4est
    displace_mesh();

    //After changing the mesh we need to rebuild things
    sim.rebuild_stokes_matrix = sim.rebuild_stokes_preconditioner = true;
    sim.computing_timer.exit_section("FreeSurface");
  }



  template <int dim>
  void Simulator<dim>::FreeSurfaceHandler::make_constraints()
  {
    if (!sim.parameters.free_surface_enabled)
      return;

    //We would like to make sure that the mesh stays conforming upon
    //redistribution, so we construct mesh_vertex_constraints, which
    //keeps track of hanging node constraints.
    mesh_vertex_constraints.clear();
    mesh_vertex_constraints.reinit(mesh_locally_relevant);

    DoFTools::make_hanging_node_constraints(free_surface_dof_handler, mesh_vertex_constraints);

    //We can safely close this now
    mesh_vertex_constraints.close();

    //Now construct the mesh displacement constraints
    mesh_displacement_constraints.clear();
    mesh_displacement_constraints.reinit(mesh_locally_relevant);

    //mesh_displacement_constraints can use the same hanging node
    //information that was used for mesh_vertex constraints.
    mesh_displacement_constraints.merge(mesh_vertex_constraints);

    //Add the vanilla periodic boundary constraints
    typedef std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int> > periodic_boundary_pairs;
    periodic_boundary_pairs pbp = sim.geometry_model->get_periodic_boundary_pairs();
    for (periodic_boundary_pairs::iterator p = pbp.begin(); p != pbp.end(); ++p)
      DoFTools::make_periodicity_constraints(free_surface_dof_handler, (*p).first.first, (*p).first.second, (*p).second, mesh_displacement_constraints);

    //Zero out the displacement for the zero-velocity boundary indicators
    for (std::set<types::boundary_id>::const_iterator p = sim.parameters.zero_velocity_boundary_indicators.begin();
         p != sim.parameters.zero_velocity_boundary_indicators.end(); ++p)
      VectorTools::interpolate_boundary_values (free_surface_dof_handler, *p,
                                                ZeroFunction<dim>(dim), mesh_displacement_constraints);

    //Zero out the displacement for the prescribed velocity
    for (std::map<types::boundary_id, std::pair<std::string, std::string> >::const_iterator p = sim.parameters.prescribed_velocity_boundary_indicators.begin();
         p != sim.parameters.prescribed_velocity_boundary_indicators.end(); ++p)
      VectorTools::interpolate_boundary_values (free_surface_dof_handler, p->first,
                                                ZeroFunction<dim>(dim), mesh_displacement_constraints);

    //make the tangential boundary indicators no displacement normal to the boundary
    VectorTools::compute_no_normal_flux_constraints (free_surface_dof_handler,
                                                     /* first_vector_component= */
                                                     0,
                                                     sim.parameters.tangential_velocity_boundary_indicators,
                                                     mesh_displacement_constraints, sim.mapping);

    //make the periodic boundary indicators no displacement normal to the boundary
    std::set< types::boundary_id > periodic_boundaries;
    for (periodic_boundary_pairs::iterator p = pbp.begin(); p != pbp.end(); ++p)
      {
        periodic_boundaries.insert((*p).first.first);
        periodic_boundaries.insert((*p).first.second);
      }
    VectorTools::compute_no_normal_flux_constraints (free_surface_dof_handler,
                                                     /* first_vector_component= */
                                                     0,
                                                     periodic_boundaries,
                                                     mesh_displacement_constraints, sim.mapping);

    //For the free surface indicators we constrain the displacement to be v.n
    LinearAlgebra::Vector boundary_normal_velocity;
    boundary_normal_velocity.reinit(mesh_locally_owned, mesh_locally_relevant, sim.mpi_communicator);
    project_normal_velocity_onto_boundary( boundary_normal_velocity );

    //now insert the relevant part of the solution into the mesh constraints
    IndexSet constrained_dofs;
    DoFTools::extract_boundary_dofs(free_surface_dof_handler, ComponentMask(dim, true),
                                    constrained_dofs, sim.parameters.free_surface_boundary_indicators);
    for ( unsigned int i = 0; i < constrained_dofs.n_elements();  ++i)
      {
        types::global_dof_index index = constrained_dofs.nth_index_in_set(i);
        if (mesh_displacement_constraints.can_store_line(index))
          if (mesh_displacement_constraints.is_constrained(index)==false)
            {
              mesh_displacement_constraints.add_line(index);
              mesh_displacement_constraints.set_inhomogeneity(index, boundary_normal_velocity[index]);
            }
      }

    mesh_displacement_constraints.close();
  }


  template <int dim>
  void Simulator<dim>::FreeSurfaceHandler::project_normal_velocity_onto_boundary(LinearAlgebra::Vector &output)
  {
    // TODO: should we use the extrapolated solution?

    //stuff for iterating over the mesh
    QGauss<dim-1> face_quadrature(free_surface_fe.degree+1);
    UpdateFlags update_flags = UpdateFlags(update_values | update_normal_vectors | update_JxW_values);
    FEFaceValues<dim> fs_fe_face_values (sim.mapping, free_surface_fe, face_quadrature, update_flags);
    FEFaceValues<dim> fe_face_values (sim.mapping, sim.finite_element, face_quadrature, update_flags);
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
    periodic_boundary_pairs pbp = sim.geometry_model->get_periodic_boundary_pairs();
    for (periodic_boundary_pairs::iterator p = pbp.begin(); p != pbp.end(); ++p)
      DoFTools::make_periodicity_constraints(free_surface_dof_handler,
                                             (*p).first.first, (*p).first.second, (*p).second, mass_matrix_constraints);

    mass_matrix_constraints.close();

    //set up the matrix
    LinearAlgebra::SparseMatrix mass_matrix;
#ifdef ASPECT_USE_PETSC
    CompressedSimpleSparsityPattern sp(mesh_locally_relevant);

#else
    TrilinosWrappers::SparsityPattern sp (mesh_locally_owned,mesh_locally_owned,
                                          sim.mpi_communicator);
#endif
    DoFTools::make_sparsity_pattern (free_surface_dof_handler, sp, mass_matrix_constraints, false,
                                     Utilities::MPI::this_mpi_process(sim.mpi_communicator));
#ifdef ASPECT_USE_PETSC
    SparsityTools::distribute_sparsity_pattern(sp,
                                               free_surface_dof_handler.n_locally_owned_dofs_per_processor(),
                                               sim.mpi_communicator, mesh_locally_relevant);

    sp.compress();
    mass_matrix.reinit (mesh_locally_owned, mesh_locally_owned, sp, sim.mpi_communicator);
#else
    sp.compress();
    mass_matrix.reinit (sp);
#endif


    FEValuesExtractors::Vector extract_vel(0);

    //make distributed vectors.
    LinearAlgebra::Vector rhs, dist_solution;
    rhs.reinit(mesh_locally_owned, sim.mpi_communicator);
    dist_solution.reinit(mesh_locally_owned, sim.mpi_communicator);

    typename DoFHandler<dim>::active_cell_iterator
    cell = sim.dof_handler.begin_active(), endc= sim.dof_handler.end();
    typename DoFHandler<dim>::active_cell_iterator
    fscell = free_surface_dof_handler.begin_active();

    for (; cell!=endc; ++cell, ++fscell)
      if (cell->at_boundary() && cell->is_locally_owned())
        for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
          if ( cell->face(face_no)->at_boundary() &&
               ((sim.parameters.free_surface_boundary_indicators.find(cell->face(face_no)->boundary_indicator())
                 != sim.parameters.free_surface_boundary_indicators.end())))
            {
              fscell->get_dof_indices (cell_dof_indices);
              fs_fe_face_values.reinit (fscell, face_no);
              fe_face_values.reinit (cell, face_no);
              fe_face_values[sim.introspection.extractors.velocities].get_function_values(sim.solution, velocity_values);

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

    //Jacobi seems to be fine here.  Other preconditioners (ILU, IC) run into troubles
    //because the matrrix is mostly empty, since we don't touch internal vertices.
    LinearAlgebra::PreconditionJacobi preconditioner_mass;
    preconditioner_mass.initialize(mass_matrix);

    SolverControl solver_control(5*rhs.size(), 1e-7*rhs.l2_norm());
    SolverCG<LinearAlgebra::Vector> cg(solver_control);
    cg.solve (mass_matrix, dist_solution, rhs, preconditioner_mass);

    mass_matrix_constraints.distribute (dist_solution);
    output = dist_solution;
  }


  template <int dim>
  void Simulator<dim>::FreeSurfaceHandler::solve_elliptic_problem()
  {
    QGauss<dim> quadrature(free_surface_fe.degree + 1);
    UpdateFlags update_flags = UpdateFlags(update_values | update_JxW_values | update_gradients);
    FEValues<dim> fe_values (sim.mapping, free_surface_fe, quadrature, update_flags);

    const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                       dofs_per_face = sim.finite_element.dofs_per_face,
                       n_q_points    = fe_values.n_quadrature_points;

    std::vector<unsigned int> cell_dof_indices (dofs_per_cell);
    std::vector<unsigned int> face_dof_indices (dofs_per_face);
    Vector<double> cell_vector (dofs_per_cell);
    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

    mesh_matrix = 0;

    //carry out the solution
    FEValuesExtractors::Vector extract_vel(0);

    LinearAlgebra::Vector rhs, poisson_solution;
    rhs.reinit(mesh_locally_owned, sim.mpi_communicator);
    poisson_solution.reinit(mesh_locally_owned, sim.mpi_communicator);

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

          mesh_displacement_constraints.distribute_local_to_global (cell_matrix, cell_vector,
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
#ifdef ASPECT_USE_PETSC
    Amg_data.symmetric_operator = false;
#else
    Amg_data.constant_modes = constant_modes;
    Amg_data.elliptic = true;
    Amg_data.higher_order_elements = false;
    Amg_data.smoother_sweeps = 2;
    Amg_data.aggregation_threshold = 0.02;
#endif
    preconditioner_stiffness.initialize(mesh_matrix);

    SolverControl solver_control(5*rhs.size(), sim.parameters.linear_stokes_solver_tolerance*rhs.l2_norm());
    SolverCG<LinearAlgebra::Vector> cg(solver_control);

    cg.solve (mesh_matrix, poisson_solution, rhs, preconditioner_stiffness);
    sim.pcout << "   Solving mesh velocity system... " << solver_control.last_step() <<" iterations."<< std::endl;

    mesh_displacement_constraints.distribute (poisson_solution);
    mesh_vertex_velocity = poisson_solution;
  }


  template <int dim>
  void Simulator<dim>::FreeSurfaceHandler::calculate_mesh_displacement()
  {
    LinearAlgebra::Vector distributed_mesh_vertices(mesh_locally_owned, sim.mpi_communicator);
    LinearAlgebra::Vector distributed_mesh_vertex_velocity(mesh_locally_owned, sim.mpi_communicator);

    distributed_mesh_vertices = mesh_vertices;
    distributed_mesh_vertex_velocity = mesh_vertex_velocity;

    //actually do the ALE thing
    distributed_mesh_vertices.sadd(1.0, sim.time_step, distributed_mesh_vertex_velocity);
    mesh_vertices = distributed_mesh_vertices;


    //Interpolate the mesh vertex velocity onto the Stokes velocity system for use in ALE corrections

    LinearAlgebra::BlockVector distributed_mesh_velocity;
    distributed_mesh_velocity.reinit(sim.introspection.index_sets.system_partitioning, sim.mpi_communicator);

    const std::vector<Point<dim> > support_points
      = sim.finite_element.base_element(sim.introspection.component_indices.velocities[0]).get_unit_support_points();

    Quadrature<dim> quad(support_points);
    UpdateFlags update_flags = UpdateFlags(update_values | update_JxW_values);
    FEValues<dim> fs_fe_values (sim.mapping, free_surface_fe, quad, update_flags);
    FEValues<dim> fe_values (sim.mapping, sim.finite_element, quad, update_flags);
    const unsigned int n_q_points = fe_values.n_quadrature_points,
                       dofs_per_cell = fe_values.dofs_per_cell;

    std::vector<unsigned int> cell_dof_indices (dofs_per_cell);
    FEValuesExtractors::Vector extract_vel(0);
    std::vector<Tensor<1,dim> > velocity_values(n_q_points);

    typename DoFHandler<dim>::active_cell_iterator
    cell = sim.dof_handler.begin_active(), endc= sim.dof_handler.end();
    typename DoFHandler<dim>::active_cell_iterator
    fscell = free_surface_dof_handler.begin_active();

    for (; cell!=endc; ++cell, ++fscell)
      if (cell->is_locally_owned())
        {
          cell->get_dof_indices (cell_dof_indices);

          fe_values.reinit (cell);
          fs_fe_values.reinit (fscell);
          fs_fe_values[extract_vel].get_function_values(mesh_vertex_velocity, velocity_values);
          for (unsigned int j=0; j<n_q_points; ++j)
            for (unsigned int dir=0; dir<dim; ++dir)
              {
                unsigned int support_point_index
                  = sim.finite_element.component_to_system_index(/*velocity component=*/ sim.introspection.component_indices.velocities[dir],
                                                                                         /*dof index within component=*/ j);
                distributed_mesh_velocity[cell_dof_indices[support_point_index]] = velocity_values[j][dir];
              }
        }

    distributed_mesh_velocity.compress(VectorOperation::insert);
    mesh_velocity = distributed_mesh_velocity;
  }


  template <int dim>
  void Simulator<dim>::FreeSurfaceHandler::setup_dofs()
  {
    if (!sim.parameters.free_surface_enabled)
      return;

    // these live in the same FE as the velocity variable:
    mesh_velocity.reinit(sim.introspection.index_sets.system_partitioning,
                         sim.introspection.index_sets.system_relevant_partitioning,
                         sim.mpi_communicator);


    free_surface_dof_handler.distribute_dofs(free_surface_fe);

    sim.pcout << "Number of free surface degrees of freedom: "
              << free_surface_dof_handler.n_dofs()
              << std::endl;

    // Renumber the DoFs hierarchical so that we get the
    // same numbering if we resume the computation. This
    // is because the numbering depends on the order the
    // cells are created.
    DoFRenumbering::hierarchical (free_surface_dof_handler);



    mesh_locally_owned = free_surface_dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs (free_surface_dof_handler,
                                             mesh_locally_relevant);

    mesh_vertices.reinit(mesh_locally_owned, mesh_locally_relevant, sim.mpi_communicator);
    mesh_vertex_velocity.reinit(mesh_locally_owned, mesh_locally_relevant, sim.mpi_communicator);

    //if we are just starting, we need to initialize mesh_vertices
    if (sim.timestep_number == 0)
      {
        LinearAlgebra::Vector distributed_mesh_vertices;
        distributed_mesh_vertices.reinit(mesh_locally_owned, sim.mpi_communicator);

        const std::vector<Point<dim> > mesh_support_points
          = free_surface_fe.base_element(0).get_unit_support_points();
        FEValues<dim> mesh_points (sim.mapping, free_surface_fe,
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

    make_constraints();

    // matrix
    {
      mesh_matrix.clear ();

      Table<2,DoFTools::Coupling> coupling (dim, dim);
      coupling.fill(DoFTools::none);

      for (unsigned int c=0; c<dim; ++c)
        coupling[c][c] = DoFTools::always;

#ifdef ASPECT_USE_PETSC
      CompressedSimpleSparsityPattern sp(mesh_locally_relevant);

#else
      TrilinosWrappers::SparsityPattern sp (mesh_locally_owned,mesh_locally_owned,
                                            sim.mpi_communicator);
#endif

      DoFTools::make_sparsity_pattern (free_surface_dof_handler,
                                       coupling, sp,
                                       mesh_displacement_constraints, false,
                                       Utilities::MPI::
                                       this_mpi_process(sim.mpi_communicator));

#ifdef ASPECT_USE_PETSC
      SparsityTools::distribute_sparsity_pattern(sp,
                                                 free_surface_dof_handler.n_locally_owned_dofs_per_processor(),
                                                 sim.mpi_communicator, mesh_locally_relevant);

      sp.compress();

      mesh_matrix.reinit (mesh_locally_owned, mesh_locally_owned, sp, sim.mpi_communicator);
#else
      sp.compress();

      mesh_matrix.reinit (sp);
#endif

    }


  }

  template <int dim>
  void Simulator<dim>::FreeSurfaceHandler::displace_mesh()
  {
    if (!sim.parameters.free_surface_enabled)
      return;

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

  template <int dim>
  void Simulator<dim>::FreeSurfaceHandler::apply_stabilization(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    FullMatrix<double> &local_matrix)
  {
    if (!sim.parameters.free_surface_enabled)
      return;

    QGauss<dim-1> quadrature(sim.parameters.stokes_velocity_degree+1);
    UpdateFlags update_flags = UpdateFlags(update_values | update_normal_vectors |
                                           update_quadrature_points | update_JxW_values);
    FEFaceValues<dim> fe_face_values (sim.mapping, sim.finite_element, quadrature, update_flags);
    const unsigned int n_face_q_points = fe_face_values.n_quadrature_points,
                       dofs_per_cell = fe_face_values.dofs_per_cell;

    const FEValuesExtractors::Vector velocities (0);

    //only apply on free surface faces
    if (cell->at_boundary() && cell->is_locally_owned())
      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
        if ( cell->face(face_no)->at_boundary() &&
             sim.parameters.free_surface_boundary_indicators.find(cell->face(face_no)->boundary_indicator())!=
             sim.parameters.free_surface_boundary_indicators.end())
          {
            fe_face_values.reinit(cell, face_no);

            //come up with the density contrast across the free surface
            std::vector<Point<dim> > quad_points = fe_face_values.get_quadrature_points();

            for (unsigned int q_point = 0; q_point < n_face_q_points; ++q_point)
              for (unsigned int i=0; i< fe_face_values.dofs_per_cell; ++i)
                for (unsigned int j=0; j< fe_face_values.dofs_per_cell; ++j)
                  {
                    //see Kaus et al 2010 for details

                    const Tensor<1,dim> gravity = sim.gravity_model->gravity_vector(quad_points[q_point]);
                    double g_norm = gravity.norm();

                    //construct the relevant vectors
                    const Tensor<1,dim> v =     fe_face_values[velocities].value(j, q_point);
                    const Tensor<1,dim> n_hat = fe_face_values.normal_vector(q_point);
                    const Tensor<1,dim> w =     fe_face_values[velocities].value(i, q_point);
                    const Tensor<1,dim> g_hat = (g_norm == 0.0 ? Tensor<1,dim>() : gravity/g_norm);

                    //TODO we should use the actual density, not the reference density here.
                    double pressure_perturbation = std::abs(sim.material_model->reference_density())*
                                                   sim.time_step*free_surface_theta*g_norm;

                    //The fictive stabilization stress is (w.g)*(v.n)
                    const double stress_value = -pressure_perturbation*
                                                (w*g_hat) * (v*n_hat)
                                                *fe_face_values.JxW(q_point);

                    local_matrix(i,j) += stress_value;
                  }
          }

  }


}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template class Simulator<dim>::FreeSurfaceHandler;

  ASPECT_INSTANTIATE(INSTANTIATE)
}
