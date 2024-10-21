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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
 */


#include <aspect/mesh_deformation/free_surface.h>
#include <aspect/simulator_signals.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/simulator/assemblers/interface.h>
#include <aspect/melt.h>
#include <aspect/simulator.h>

#include <deal.II/dofs/dof_tools.h>


namespace aspect
{
  namespace MeshDeformation
  {
    template <int dim>
    void
    FreeSurface<dim>::initialize ()
    {
      // Pressure normalization doesn't really make sense with a free surface, and if we do
      // use it, we can run into problems with geometry_model->depth().
      AssertThrow ( this->get_parameters().pressure_normalization == "no",
                    ExcMessage("The free surface scheme can only be used with no pressure normalization") );

      // Check that we do not use the free surface on a boundary that has zero slip,
      // free slip or prescribed velocity boundary conditions on it.

      // Get the zero velocity boundary indicators
      std::set<types::boundary_id> velocity_boundary_indicators = this->get_boundary_velocity_manager().get_zero_boundary_velocity_indicators();

      // Get the tangential velocity boundary indicators
      const auto &tangential_boundary_indicators = this->get_boundary_velocity_manager().get_tangential_boundary_velocity_indicators();
      velocity_boundary_indicators.insert(tangential_boundary_indicators.begin(),
                                          tangential_boundary_indicators.end());

      // Get the active velocity boundary indicators
      const auto &prescribed_boundary_indicators = this->get_boundary_velocity_manager().get_prescribed_boundary_velocity_indicators();
      velocity_boundary_indicators.insert(prescribed_boundary_indicators.begin(),
                                          prescribed_boundary_indicators.end());

      // Get the mesh deformation boundary indicators
      for (const auto &p : this->get_mesh_deformation_boundary_indicators())
        AssertThrow(velocity_boundary_indicators.find(p) == velocity_boundary_indicators.end(),
                    ExcMessage("The free surface mesh deformation plugin cannot be used with the current velocity boundary conditions"));
    }


    template <int dim>
    void FreeSurface<dim>::project_velocity_onto_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                          const IndexSet &mesh_locally_owned,
                                                          const IndexSet &mesh_locally_relevant,
                                                          LinearAlgebra::Vector &output) const
    {
      // TODO: should we use the extrapolated solution?

      // stuff for iterating over the mesh
      const QGauss<dim-1> face_quadrature(mesh_deformation_dof_handler.get_fe().degree+1);
      UpdateFlags update_flags = UpdateFlags(update_values | update_quadrature_points
                                             | update_normal_vectors | update_JxW_values);
      FEFaceValues<dim> fs_fe_face_values (this->get_mapping(), mesh_deformation_dof_handler.get_fe(), face_quadrature, update_flags);
      FEFaceValues<dim> fe_face_values (this->get_mapping(), this->get_fe(), face_quadrature, update_flags);
      const unsigned int n_face_q_points = fe_face_values.n_quadrature_points,
                         dofs_per_cell = fs_fe_face_values.dofs_per_cell;

      // stuff for assembling system
      std::vector<types::global_dof_index> cell_dof_indices (dofs_per_cell);
      Vector<double> cell_vector (dofs_per_cell);
      FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

      // stuff for getting the velocity values
      std::vector<Tensor<1,dim>> velocity_values(n_face_q_points);

      // set up constraints
      AffineConstraints<double> mass_matrix_constraints(
#if DEAL_II_VERSION_GTE(9,6,0)
        mesh_deformation_dof_handler.locally_owned_dofs(),
#endif
        mesh_locally_relevant);
      DoFTools::make_hanging_node_constraints(mesh_deformation_dof_handler, mass_matrix_constraints);

      using periodic_boundary_pairs = std::set<std::pair<std::pair<types::boundary_id, types::boundary_id>, unsigned int>>;
      periodic_boundary_pairs pbp = this->get_geometry_model().get_periodic_boundary_pairs();
      for (const auto &p : pbp)
        DoFTools::make_periodicity_constraints(mesh_deformation_dof_handler,
                                               p.first.first, p.first.second, p.second, mass_matrix_constraints);

      mass_matrix_constraints.close();

      // set up the matrix
      LinearAlgebra::SparseMatrix mass_matrix;
      TrilinosWrappers::SparsityPattern sp (mesh_locally_owned,
                                            mesh_locally_owned,
                                            mesh_locally_relevant,
                                            this->get_mpi_communicator());
      DoFTools::make_sparsity_pattern (mesh_deformation_dof_handler, sp, mass_matrix_constraints, false,
                                       Utilities::MPI::this_mpi_process(this->get_mpi_communicator()));
      sp.compress();
      mass_matrix.reinit (sp);

      FEValuesExtractors::Vector extract_vel(0);

      // make distributed vectors.
      LinearAlgebra::Vector rhs, dist_solution;
      rhs.reinit(mesh_locally_owned, this->get_mpi_communicator());
      dist_solution.reinit(mesh_locally_owned, this->get_mpi_communicator());

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(), endc= this->get_dof_handler().end();
      typename DoFHandler<dim>::active_cell_iterator
      fscell = mesh_deformation_dof_handler.begin_active();

      // Get the boundary indicators of those boundaries with
      // a free surface.
      const std::set<types::boundary_id> tmp_free_surface_boundary_indicators = this->get_mesh_deformation_handler().get_free_surface_boundary_indicators();

      for (; cell!=endc; ++cell, ++fscell)
        if (cell->at_boundary() && cell->is_locally_owned())
          for (const unsigned int face_no : cell->face_indices())
            if (cell->face(face_no)->at_boundary())
              {
                const types::boundary_id boundary_indicator
                  = cell->face(face_no)->boundary_id();

                // Only project onto the free surface boundary/boundaries.
                if (tmp_free_surface_boundary_indicators.find(boundary_indicator) == tmp_free_surface_boundary_indicators.end())
                  continue;

                fscell->get_dof_indices (cell_dof_indices);
                fs_fe_face_values.reinit (fscell, face_no);
                fe_face_values.reinit (cell, face_no);
                fe_face_values[this->introspection().extractors.velocities].get_function_values(this->get_solution(), velocity_values);

                cell_vector = 0;
                cell_matrix = 0;
                for (unsigned int q=0; q<n_face_q_points; ++q)
                  {
                    // Select the direction onto which to project the velocity solution
                    Tensor<1,dim> direction;
                    if ( advection_direction == SurfaceAdvection::normal ) // project onto normal vector
                      direction = fs_fe_face_values.normal_vector(q);
                    else if ( advection_direction == SurfaceAdvection::vertical ) // project onto local gravity
                      direction = this->get_gravity_model().gravity_vector(fs_fe_face_values.quadrature_point(q));
                    else
                      AssertThrow(false, ExcInternalError());

                    direction *= ( direction.norm() > 0.0 ? 1./direction.norm() : 0.0 );

                    const double JxW = fs_fe_face_values.JxW(q);

                    small_vector<Tensor<1,dim>> phi_u(dofs_per_cell);
                    for (unsigned int i=0; i<dofs_per_cell; ++i)
                      phi_u[i] = fs_fe_face_values[extract_vel].value(i,q);

                    for (unsigned int i=0; i<dofs_per_cell; ++i)
                      {
                        for (unsigned int j=0; j<dofs_per_cell; ++j)
                          {
                            cell_matrix(i,j) += (phi_u[j] *
                                                 phi_u[i]) *
                                                JxW;
                          }

                        cell_vector(i) += (phi_u[i] * direction)
                                          * (velocity_values[q] * direction)
                                          * JxW;
                      }
                  }

                mass_matrix_constraints.distribute_local_to_global (cell_matrix, cell_vector,
                                                                    cell_dof_indices, mass_matrix, rhs, false);
              }

      rhs.compress (VectorOperation::add);
      mass_matrix.compress(VectorOperation::add);

      // Jacobi seems to be fine here.  Other preconditioners (ILU, IC) run into troubles
      // because the matrix is mostly empty, since we don't touch internal vertices.
      LinearAlgebra::PreconditionJacobi preconditioner_mass;
      preconditioner_mass.initialize(mass_matrix);

      SolverControl solver_control(5*rhs.size(), this->get_parameters().linear_stokes_solver_tolerance*rhs.l2_norm());
      SolverCG<LinearAlgebra::Vector> cg(solver_control);
      cg.solve (mass_matrix, dist_solution, rhs, preconditioner_mass);

      mass_matrix_constraints.distribute (dist_solution);
      output = dist_solution;
    }



    /**
     * A function that creates constraints for the velocity of certain mesh
     * vertices (e.g. the surface vertices) for a specific boundary.
     * The calling class will respect
     * these constraints when computing the new vertex positions.
     */
    template <int dim>
    void
    FreeSurface<dim>::compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                               AffineConstraints<double> &mesh_velocity_constraints,
                                                               const std::set<types::boundary_id> &boundary_id) const
    {
      // For the free surface indicators we constrain the displacement to be v.n
      LinearAlgebra::Vector boundary_velocity;

      const IndexSet &mesh_locally_owned = mesh_deformation_dof_handler.locally_owned_dofs();
      const IndexSet mesh_locally_relevant = DoFTools::extract_locally_relevant_dofs (mesh_deformation_dof_handler);
      boundary_velocity.reinit(mesh_locally_owned, mesh_locally_relevant,
                               this->get_mpi_communicator());
      project_velocity_onto_boundary(mesh_deformation_dof_handler, mesh_locally_owned,
                                     mesh_locally_relevant, boundary_velocity);

      // now insert the relevant part of the solution into the mesh constraints
      const IndexSet constrained_dofs =
        DoFTools::extract_boundary_dofs(mesh_deformation_dof_handler,
                                        ComponentMask(dim, true),
                                        boundary_id);

      for (unsigned int i = 0; i < constrained_dofs.n_elements();  ++i)
        {
          types::global_dof_index index = constrained_dofs.nth_index_in_set(i);
          if (mesh_velocity_constraints.can_store_line(index))
            if (mesh_velocity_constraints.is_constrained(index)==false)
              {
#if DEAL_II_VERSION_GTE(9,6,0)
                mesh_velocity_constraints.add_constraint(index,
                                                         {},
                                                         boundary_velocity[index]);
#else
                mesh_velocity_constraints.add_line(index);
                mesh_velocity_constraints.set_inhomogeneity(index, boundary_velocity[index]);
#endif
              }
        }
    }



    template <int dim>
    bool
    FreeSurface<dim>::
    needs_surface_stabilization () const
    {
      return true;
    }



    template <int dim>
    void FreeSurface<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection ("Free surface");
        {
          prm.declare_entry("Surface velocity projection", "normal",
                            Patterns::Selection("normal|vertical"),
                            "After each time step the free surface must be "
                            "advected in the direction of the velocity field. "
                            "Mass conservation requires that the mesh velocity "
                            "is in the normal direction of the surface. However, "
                            "for steep topography or large curvature, advection "
                            "in the normal direction can become ill-conditioned, "
                            "and instabilities in the mesh can form. Projection "
                            "of the mesh velocity onto the local vertical direction "
                            "can preserve the mesh quality better, but at the "
                            "cost of slightly poorer mass conservation of the "
                            "domain.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }



    template <int dim>
    void FreeSurface<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection ("Free surface");
        {
          std::string advection_dir = prm.get("Surface velocity projection");

          if ( advection_dir == "normal")
            advection_direction = SurfaceAdvection::normal;
          else if ( advection_dir == "vertical")
            advection_direction = SurfaceAdvection::vertical;
          else
            AssertThrow(false, ExcMessage("The surface velocity projection must be ``normal'' or ``vertical''."));
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace MeshDeformation
  {
    ASPECT_REGISTER_MESH_DEFORMATION_MODEL(FreeSurface,
                                           "free surface",
                                           "A plugin that computes the deformation of surface "
                                           "vertices according to the solution of the flow problem. "
                                           "In particular this means if the surface of the domain is "
                                           "left open to flow, this flow will carry the mesh with it. "
                                           "The implementation was described in \\cite{rose_freesurface}, "
                                           "with the stabilization of the free surface originally described "
                                           "in \\cite{kaus:etal:2010}.")
  }
}
