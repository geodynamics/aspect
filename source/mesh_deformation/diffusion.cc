/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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


#include <aspect/mesh_deformation/diffusion.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/two_merged_boxes.h>
#include <aspect/simulator.h>
#include <aspect/geometry_model/initial_topography_model/zero_topography.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <boost/lexical_cast.hpp>

namespace aspect
{

  namespace MeshDeformation
  {
    template <int dim>
    Diffusion<dim>::Diffusion()
      :
      diffusivity(0),
      start_time(std::numeric_limits<double>::quiet_NaN()),
      current_time(0),
      last_diffusion_time(std::numeric_limits<double>::quiet_NaN()),
      time_between_diffusion(std::numeric_limits<double>::max()),
      start_timestep(0),
      last_diffusion_timestep (1),
      timesteps_between_diffusion (1),
      apply_diffusion(false)
    {}



    template <int dim>
    void
    Diffusion<dim>::initialize ()
    {
      AssertThrow(Plugins::plugin_type_matches<GeometryModel::Box<dim> >(this->get_geometry_model()) ||
                  Plugins::plugin_type_matches<GeometryModel::TwoMergedBoxes<dim> >(this->get_geometry_model()),
                  ExcMessage("The surface diffusion mesh deformation plugin only works for Box geometries. "));

    }



    template <int dim>
    void
    Diffusion<dim>::update ()
    {
      // Initialize the start time and timestep
      if (std::isnan(start_time))
        {
          start_time = (this->get_time());
          start_timestep = this->get_timestep_number();
        }

      // Set the current time and timestep
      current_time = (this->get_time());
      const unsigned int current_timestep_number = this->get_timestep_number();

      // Determine whether we need to apply diffusion based
      // on the time or timestep interval between applications.
      if (current_timestep_number != 0)
        {
          if ((std::isnan(last_diffusion_time) && current_time >= start_time + time_between_diffusion)
              || (std::isnan(last_diffusion_time) && current_timestep_number >= start_timestep + timesteps_between_diffusion))
            apply_diffusion = true;
          else if ((current_time >= last_diffusion_time + time_between_diffusion)
                   || (current_timestep_number >= last_diffusion_timestep + timesteps_between_diffusion))
            apply_diffusion = true;
          else
            apply_diffusion = false;
        }

      // If diffusion is applied, set the current time(step) as
      // the last time(step) of application.
      if (apply_diffusion)
        {
          last_diffusion_time = current_time;
          last_diffusion_timestep = current_timestep_number;
        }

    }



    template <int dim>
    void Diffusion<dim>::diffuse_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                          const IndexSet &mesh_locally_owned,
                                          const IndexSet &mesh_locally_relevant,
                                          LinearAlgebra::Vector &output,
                                          const std::set<types::boundary_id> boundary_ids) const
    {
      // Check that the current timestep does not exceed the diffusion timestep
      compute_time_step(mesh_deformation_dof_handler, boundary_ids);

      // Set up constraints
      ConstraintMatrix mass_matrix_constraints(mesh_locally_relevant);
      DoFTools::make_hanging_node_constraints(mesh_deformation_dof_handler, mass_matrix_constraints);

      typedef std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int> > periodic_boundary_pairs;
      periodic_boundary_pairs pbp = this->get_geometry_model().get_periodic_boundary_pairs();
      for (periodic_boundary_pairs::iterator p = pbp.begin(); p != pbp.end(); ++p)
        DoFTools::make_periodicity_constraints(mesh_deformation_dof_handler,
                                               (*p).first.first, (*p).first.second, (*p).second, mass_matrix_constraints);

      // The list of boundary indicators fow which we need to set
      // zero mesh velocities, which means the zero/prescribed Stokes velocity
      // boundaries minus those that are listed as Additional
      // tangential boundary indicators.
      std::set<types::boundary_id> x_zero_boundary_indicators = zero_boundary_velocity_indicators;
      x_zero_boundary_indicators.insert(prescribed_boundary_velocity_indicators.begin(), prescribed_boundary_velocity_indicators.end());
      for (std::set<types::boundary_id>::const_iterator p = x_zero_boundary_indicators.begin();
           p != x_zero_boundary_indicators.end(); ++p)
        if (boundary_ids.find(*p) == boundary_ids.end())
          if (additional_tangential_mesh_boundary_indicators.find(*p) == additional_tangential_mesh_boundary_indicators.end())
            {
              VectorTools::interpolate_boundary_values (mesh_deformation_dof_handler, *p,
                                                        ZeroFunction<dim>(dim), mass_matrix_constraints);
            }

      // The list of boundary indicators for which we need to set
      // no_normal_flux_constraints, which means all
      // minus the diffusion mesh deformation boundary indicators
      // and minus the zero/prescribed Stokes velocity boundary indicators.
      std::set<types::boundary_id> x_no_flux_boundary_indicators = tangential_boundary_velocity_indicators;
      x_no_flux_boundary_indicators.insert(additional_tangential_mesh_boundary_indicators.begin(),additional_tangential_mesh_boundary_indicators.end());
      for (std::set<types::boundary_id>::const_iterator p = x_no_flux_boundary_indicators.begin();
           p != x_no_flux_boundary_indicators.end(); ++p)
        if (boundary_ids.find(*p) != boundary_ids.end())
          {
            x_no_flux_boundary_indicators.erase(*p);
          }

      // Make the no flux boundary constraints
      VectorTools::compute_no_normal_flux_constraints (mesh_deformation_dof_handler,
                                                       /* first_vector_component= */
                                                       0,
                                                       x_no_flux_boundary_indicators,
                                                       mass_matrix_constraints, this->get_mapping());

      mass_matrix_constraints.close();

      // Set up the system to solve
      LinearAlgebra::SparseMatrix mass_matrix;

      // Sparsity of the matrix
#ifdef ASPECT_USE_PETSC
      LinearAlgebra::DynamicSparsityPattern sp(mesh_locally_relevant);

#else
      TrilinosWrappers::SparsityPattern sp (mesh_locally_owned,
                                            mesh_locally_owned,
                                            mesh_locally_relevant,
                                            this->get_mpi_communicator());
#endif
      DoFTools::make_sparsity_pattern (mesh_deformation_dof_handler, sp, mass_matrix_constraints, false,
                                       Utilities::MPI::this_mpi_process(this->get_mpi_communicator()));

#ifdef ASPECT_USE_PETSC
      SparsityTools::distribute_sparsity_pattern(sp,
                                                 mesh_deformation_dof_handler.n_locally_owned_dofs_per_processor(),
                                                 this->get_mpi_communicator(), mesh_locally_relevant);

      sp.compress();
      mass_matrix.reinit (mesh_locally_owned, mesh_locally_owned, sp, this->get_mpi_communicator());
#else
      sp.compress();
      mass_matrix.reinit (sp);
#endif

      LinearAlgebra::Vector system_rhs, solution;
      system_rhs.reinit(mesh_locally_owned, this->get_mpi_communicator());
      solution.reinit(mesh_locally_owned, this->get_mpi_communicator());

      // Initialize Gauss-Legendre quadrature for degree+1 quadrature points of the surface faces
      QGauss<dim-1> face_quadrature(mesh_deformation_dof_handler.get_fe().degree+1);
      // Update shape function values and gradients, the quadrature points and the Jacobian x quadrature weights.
      UpdateFlags update_flags = UpdateFlags(update_values | update_gradients | update_quadrature_points | update_normal_vectors | update_JxW_values);
      // We want to extract the displacement at the free surface faces of the mesh deformation element.
      FEFaceValues<dim> fs_fe_face_values (this->get_mapping(), mesh_deformation_dof_handler.get_fe(), face_quadrature, update_flags);
      // and to solve on the whole mesh deformation mesh
      // The number of quadrature points on a mesh deformation surface face
      const unsigned int n_fs_face_q_points = fs_fe_face_values.n_quadrature_points;

      // What we need to build our system on the mesh deformation element

      // The nr of shape functions per mesh deformation element
      const unsigned int dofs_per_cell = mesh_deformation_dof_handler.get_fe().dofs_per_cell;

      // Map of local to global cell dof indices
      std::vector<types::global_dof_index> cell_dof_indices (dofs_per_cell);

      // The local rhs vector
      Vector<double> cell_vector (dofs_per_cell);
      // The local matrix
      FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

      // Vector for getting the local dim displacement values
      std::vector<Tensor<1, dim> > displacement_values(n_fs_face_q_points);

      // Vector for getting the local dim initial topography values
      std::vector<Tensor<1, dim> > initial_topography_values(n_fs_face_q_points);

      // The global displacements on the MeshDeformation FE
      LinearAlgebra::Vector displacements = this->get_mesh_deformation_handler().get_mesh_displacements();

      // The global initial topography on the MeshDeformation FE
      // TODO find another way to get to the initial topography?
      LinearAlgebra::Vector initial_topography = this->get_mesh_deformation_handler().get_initial_topography();

      // Do nothing at time zero
      if (this->get_timestep_number() < 1)
        return;

      // An extractor for the dim-valued displacement vectors
      // Later on we will compute the gravity-parallel displacement
      FEValuesExtractors::Vector extract_vertical_displacements(0);

      // An extractor for the dim-valued initial topography vectors
      // Later on we will compute the gravity-parallel displacement
      FEValuesExtractors::Vector extract_initial_topography(0);

      // Cell iterator over the MeshDeformation FE
      typename DoFHandler<dim>::active_cell_iterator
      fscell = mesh_deformation_dof_handler.begin_active(),
      fsendc= mesh_deformation_dof_handler.end();

      // Iterate over all cells to find those at the mesh deformation boundary
      for (; fscell!=fsendc; ++fscell)
        if (fscell->at_boundary() && fscell->is_locally_owned())
          for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
            if (fscell->face(face_no)->at_boundary())
              {
                // Boundary indicator of current cell face
                const types::boundary_id boundary_indicator
                  = fscell->face(face_no)->boundary_id();

                // Only apply diffusion to the requested boundaries
                if (boundary_ids.find(boundary_indicator) == boundary_ids.end())
                  continue;

                // Recompute values, gradients, etc on the faces
                fs_fe_face_values.reinit (fscell, face_no);

                // Get the global numbers of the local DoFs of the mesh deformation cell
                fscell->get_dof_indices (cell_dof_indices);

                // Extract the displacement values
                fs_fe_face_values[extract_vertical_displacements].get_function_values (displacements, displacement_values);

                // Extract the initial topography values
                fs_fe_face_values[extract_initial_topography].get_function_values (initial_topography, initial_topography_values);

                // Reset local rhs and matrix
                cell_vector = 0;
                cell_matrix = 0;

                // Loop over the quadrature points of the current face
                for (unsigned int point=0; point<n_fs_face_q_points; ++point)
                  {
                    // Get the gravity vector to compute the outward direction of displacement
                    Tensor<1,dim> direction = -(this->get_gravity_model().gravity_vector(fs_fe_face_values.quadrature_point(point)));
                    // Normalize direction vector
                    if (direction.norm() > 0.0)
                      direction *= 1./direction.norm();
                    // TODO this is only correct for box geometries
                    else
                      direction[dim-1] = 1.;

                    // Compute the total displacement in the gravity direction,
                    // i.e. the initial topography + any additional mesh displacement.
                    const double displacement = (displacement_values[point] * direction) + (initial_topography_values[point] * direction);

                    // To project onto the tangent space of the surface,
                    // we define the projection P:= I- n x n,
                    // with I the unit tensor and n the unit normal to the surface.
                    // The surface gradient then is P times the usual gradient of the shape functions.
                    const Tensor<2, dim, double> projection = unit_symmetric_tensor<dim>() -
                                                              outer_product(fs_fe_face_values.normal_vector(point), fs_fe_face_values.normal_vector(point));

                    // Loop over the shape functions
                    for (unsigned int i=0; i<dofs_per_cell; ++i)
                      {
                        // Make sure we only assemble for the y-component
                        // TODO this is only correct for box geometries
                        if (mesh_deformation_dof_handler.get_fe().system_to_component_index(i).first != dim-1)
                          continue;

                        // Assemble the RHS
                        // RHS = M*H_old
                        cell_vector(i) += displacement *
                                          fs_fe_face_values.shape_value (i, point) *
                                          fs_fe_face_values.JxW(point);


                        for (unsigned int j=0; j<dofs_per_cell; ++j)
                          {
                            // Make sure we only assemble for the y-component
                            // TODO this is only correct for box geometries
                            if (mesh_deformation_dof_handler.get_fe().system_to_component_index(j).first != dim-1)
                              continue;
                            // Assemble the matrix, for backward first order time discretization:
                            // Matrix := (M+dt*K) = (M+dt*B^T*kappa*B)
                            cell_matrix(i,j) +=
                              (
                                this->get_timestep() * diffusivity *
                                projection * fs_fe_face_values.shape_grad(i, point) * projection * fs_fe_face_values.shape_grad(j,point) +
                                fs_fe_face_values.shape_value (i, point) * fs_fe_face_values.shape_value (j, point)
                              )
                              * fs_fe_face_values.JxW(point);
                          }
                      }

                  }

                mass_matrix_constraints.distribute_local_to_global (cell_matrix, cell_vector,
                                                                    cell_dof_indices, mass_matrix, system_rhs, false);
              }

      system_rhs.compress (VectorOperation::add);
      mass_matrix.compress(VectorOperation::add);

      // Jacobi seems to be fine here.  Other preconditioners (ILU, IC) run into trouble
      // because the matrix is mostly empty, since we don't touch internal vertices.
      LinearAlgebra::PreconditionJacobi preconditioner_mass;
      LinearAlgebra::PreconditionJacobi::AdditionalData preconditioner_control(1,1e-16,1);
      preconditioner_mass.initialize(mass_matrix, preconditioner_control);

      this->get_pcout() << "   Solving mesh surface diffusion" << std::endl;
      SolverControl solver_control(50*system_rhs.size(), this->get_parameters().linear_stokes_solver_tolerance*system_rhs.l2_norm());
      SolverCG<LinearAlgebra::Vector> cg(solver_control);
      cg.solve (mass_matrix, solution, system_rhs, preconditioner_mass);

      // Distribute constraints on mass matrix
      mass_matrix_constraints.distribute (solution);

      // The solution contains the new displacements, but we need to return a velocity.
      // Therefore, we compute v=d_displacement/d_t.
      // d_displacement are the new mesh node locations
      // minus the old locations, which are initial_topography + displacements.
      LinearAlgebra::Vector d_displacement(mesh_locally_owned, mesh_locally_relevant, this->get_mpi_communicator());
      d_displacement = solution;
      d_displacement -= initial_topography;
      d_displacement -= displacements;

      // The velocity
      if (this->get_timestep() > 0.)
        d_displacement /= this->get_timestep();
      else
        d_displacement = 0.;

      output = d_displacement;
    }


    template <int dim>
    void Diffusion<dim>::compute_time_step (const DoFHandler<dim> &mesh_deformation_dof_handler,
                                            const std::set<types::boundary_id> boundary_ids) const
    {
      // Initialize Gauss-Legendre quadrature for degree+1 quadrature points of the surface faces
      const QGauss<dim-1> face_quadrature(mesh_deformation_dof_handler.get_fe().degree+1);

      // TODO Do we need to update anything for the vertex distance?
      FEFaceValues<dim> fs_fe_face_values (this->get_mapping(), mesh_deformation_dof_handler.get_fe(), face_quadrature, update_default);

      double min_local_conduction_timestep = std::numeric_limits<double>::max();

      for (const auto &fscell : mesh_deformation_dof_handler.active_cell_iterators())
        if (fscell->at_boundary() && fscell->is_locally_owned())
          for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
            if (fscell->face(face_no)->at_boundary())
              {
                // Get the boundary indicator of current cell face
                const types::boundary_id boundary_indicator
                  = fscell->face(face_no)->boundary_id();

                // Only consider the requested boundaries
                if (boundary_ids.find(boundary_indicator) == boundary_ids.end())
                  continue;

                // Reninitalize update flags for current cell face
                fs_fe_face_values.reinit (fscell, face_no);

                // Calculate the corresponding conduction timestep, if applicable
                if (diffusivity > 0.)
                  {
                    min_local_conduction_timestep = std::min(min_local_conduction_timestep,
                                                             this->get_parameters().CFL_number*pow(fscell->face(face_no)->minimum_vertex_distance(),2.)
                                                             / diffusivity);

                  }
              }

      // Get the global minimum timestep
      const double min_conduction_timestep = - Utilities::MPI::max (-min_local_conduction_timestep, this->get_mpi_communicator());

      AssertThrow (min_conduction_timestep > 0.,
                   ExcMessage("The time step length for diffusion of the surface needs to be positive, "
                              "but the computed step length was: " + std::to_string(min_conduction_timestep) + ". "
                              "Please check for non-positive diffusivity."));

      double conduction_timestep = min_conduction_timestep;
      if (this->convert_output_to_years())
        conduction_timestep /= year_in_seconds;

      AssertThrow (this->get_timestep() <= min_conduction_timestep,
                   ExcMessage("The numerical timestep is too large for diffusion of the surface. Although the "
                              "diffusion scheme is stable, note that the error increases linearly with the timestep. "
                              "The diffusion timestep is: " + std::to_string(conduction_timestep) + ". "));
    }


    /**
     * A function that creates constraints for the velocity of certain mesh
     * vertices (e.g. the surface vertices) for a specific boundary.
     * The calling class will respect
     * these constraints when computing the new vertex positions.
     */
    template <int dim>
    void
    Diffusion<dim>::compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                             ConstraintMatrix &mesh_velocity_constraints,
                                                             const std::set<types::boundary_id> &boundary_id) const
    {
      if (!apply_diffusion)
        return;

      LinearAlgebra::Vector boundary_velocity;

      const IndexSet &mesh_locally_owned = mesh_deformation_dof_handler.locally_owned_dofs();
      IndexSet mesh_locally_relevant;
      DoFTools::extract_locally_relevant_dofs (mesh_deformation_dof_handler,
                                               mesh_locally_relevant);
      boundary_velocity.reinit(mesh_locally_owned, mesh_locally_relevant,
                               this->get_mpi_communicator());

      // Determine the mesh velocity at the surface based on diffusion of
      // the topography
      diffuse_boundary(mesh_deformation_dof_handler, mesh_locally_owned,
                       mesh_locally_relevant, boundary_velocity, boundary_id);

      // now insert the relevant part of the solution into the mesh constraints
      IndexSet constrained_dofs;
      DoFTools::extract_boundary_dofs(mesh_deformation_dof_handler,
                                      ComponentMask(dim, true),
                                      constrained_dofs,
                                      boundary_id);

      for (unsigned int i = 0; i < constrained_dofs.n_elements();  ++i)
        {
          types::global_dof_index index = constrained_dofs.nth_index_in_set(i);
          if (mesh_velocity_constraints.can_store_line(index))
            if (mesh_velocity_constraints.is_constrained(index)==false)
              {
                mesh_velocity_constraints.add_line(index);
                mesh_velocity_constraints.set_inhomogeneity(index, boundary_velocity[index]);
              }
        }
    }


    template <int dim>
    void Diffusion<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh deformation");
      {
        prm.enter_subsection("Diffusion");
        {
          prm.declare_entry("Hillslope transport coefficient", "0.5",
                            Patterns::Double(0),
                            "The hillslope transport coefficient used to "
                            "diffuse the free surface, either as a  "
                            "stabilization step or a to mimic erosional "
                            "and depositional processes. Units: m2/s. ");
          prm.declare_entry("Time between diffusion", boost::lexical_cast<std::string>(std::numeric_limits<double>::max()),
                            Patterns::Double(0,std::numeric_limits<double>::max()),
                            "The time between each application of diffusion. "
                            "Units: years if the "
                            "'Use years in output instead of seconds' parameter is set; "
                            "seconds otherwise.");
          prm.declare_entry("Time steps between diffusion", "1",
                            Patterns::Integer(0,std::numeric_limits<int>::max()),
                            "The maximum number of time steps between each application of "
                            "diffusion.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void Diffusion<dim>::parse_parameters(ParameterHandler &prm)
    {
      // The list of tangential Stokes velocity boundary indicators.
      tangential_boundary_velocity_indicators = this->get_boundary_velocity_manager().get_tangential_boundary_velocity_indicators();
      // The list of zero Stokes velocity boundary indicators.
      zero_boundary_velocity_indicators = this->get_boundary_velocity_manager().get_zero_boundary_velocity_indicators();
      // The list of prescribed Stokes velocity boundary indicators.
      const std::map<types::boundary_id, std::pair<std::string,std::vector<std::string> > > active_boundary_velocity_indicators =
        this->get_boundary_velocity_manager().get_active_boundary_velocity_names();
      for (std::map<types::boundary_id, std::pair<std::string, std::vector<std::string> > >::const_iterator p = active_boundary_velocity_indicators.begin();
           p != active_boundary_velocity_indicators.end(); ++p)
        prescribed_boundary_velocity_indicators.insert(p->first);

      prm.enter_subsection ("Mesh deformation");
      {
        // Create the list of tangential mesh movement boundary indicators.
        try
          {
            const std::vector<types::boundary_id> x_additional_tangential_mesh_boundary_indicators
              = this->get_geometry_model().translate_symbolic_boundary_names_to_ids(Utilities::split_string_list
                                                                                    (prm.get ("Additional tangential mesh velocity boundary indicators")));

            additional_tangential_mesh_boundary_indicators.insert(x_additional_tangential_mesh_boundary_indicators.begin(),
                                                                  x_additional_tangential_mesh_boundary_indicators.end());
          }
        catch (const std::string &error)
          {
            AssertThrow (false, ExcMessage ("While parsing the entry <Mesh deformation/Additional tangential "
                                            "mesh velocity boundary indicators>, there was an error. Specifically, "
                                            "the conversion function complained as follows: "
                                            + error));
          }

        prm.enter_subsection ("Diffusion");
        {
          diffusivity                 = prm.get_double("Hillslope transport coefficient");
          time_between_diffusion      = prm.get_double("Time between diffusion");
          timesteps_between_diffusion = prm.get_integer("Time steps between diffusion");
          if (this->convert_output_to_years())
            time_between_diffusion *= year_in_seconds;
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
    ASPECT_REGISTER_MESH_DEFORMATION_MODEL(Diffusion,
                                           "diffusion",
                                           "A plugin that computes the deformation of surface "
                                           "vertices according to the solution of the hillslope diffusion problem. "
                                           "Diffusion can be applied every timestep, mimicking surface processing, "
                                           "or at a user-defined time or timestep interval to purely smooth the surface "
                                           "topography to avoid too great distortion of mesh elements when a free "
                                           "surface is used.")
  }
}
