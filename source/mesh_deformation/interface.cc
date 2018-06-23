/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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


#include <aspect/mesh_deformation/interface.h>
#include <aspect/simulator.h>
#include <aspect/global.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1_eulerian.h>

#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/numerics/vector_tools.h>



namespace aspect
{
  namespace MeshDeformation
  {
    template <int dim>
    Interface<dim>::~Interface ()
    {}

    template <int dim>
    void
    Interface<dim>::initialize ()
    {}

    template <int dim>
    void
    Interface<dim>::update ()
    {}

    template <int dim>
    void
    Interface<dim>::
    declare_parameters (ParameterHandler &)
    {}


    template <int dim>
    void
    Interface<dim>::parse_parameters (ParameterHandler &)
    {}


    template <int dim>
    FreeSurfaceHandler<dim>::FreeSurfaceHandler (Simulator<dim> &simulator)
      : sim(simulator),  // reference to the simulator that owns the FreeSurfaceHandler
        free_surface_fe (FE_Q<dim>(1),dim), // Q1 elements which describe the mesh geometry
        free_surface_dof_handler (sim.triangulation)
    {}

    template <int dim>
    FreeSurfaceHandler<dim>::~FreeSurfaceHandler ()
    {
      // Free the Simulator's mapping object, otherwise
      // when the FreeSurfaceHandler gets destroyed,
      // the mapping's reference to the mesh displacement
      // vector will be invalid.
      sim.mapping.reset();
    }



    namespace
    {
      std_cxx11::tuple
      <void *,
      void *,
      aspect::internal::Plugins::PluginList<Interface<2> >,
      aspect::internal::Plugins::PluginList<Interface<3> > > registered_plugins;
    }



    template <int dim>
    void
    FreeSurfaceHandler<dim>::register_mesh_deformation (const std::string &name,
                                                        const std::string &description,
                                                        void (*declare_parameters_function) (ParameterHandler &),
                                                        Interface<dim> *(*factory_function) ())
    {
      std_cxx11::get<dim>(registered_plugins).register_plugin (name,
                                                               description,
                                                               declare_parameters_function,
                                                               factory_function);
    }




    template <int dim>
    void
    FreeSurfaceHandler<dim>::update ()
    {
      for (typename std::vector<std_cxx11::shared_ptr<Interface<dim> > >::iterator
           model = mesh_deformation_objects.begin(); model != mesh_deformation_objects.end(); ++model)
        (*model)->update();
    }



    template <int dim>
    void FreeSurfaceHandler<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Free surface");
      {
        const std::string pattern_of_names
          = std_cxx11::get<dim>(registered_plugins).get_pattern_of_names ();

        prm.declare_entry("List of model names",
                          "free surface",
                          Patterns::MultipleSelection(pattern_of_names),
                          "A comma-separated list of mesh deformation models that "
                          "will be used to deform the mesh over the model runtime. "
                          "These plugins are applied in the order given. \n\n"
                          "The following mesh deformation models are available:\n\n"
                          +
                          std_cxx11::get<dim>(registered_plugins).get_description_string());

        prm.declare_entry ("Additional tangential mesh velocity boundary indicators", "",
                           Patterns::List (Patterns::Anything()),
                           "A comma separated list of names denoting those boundaries "
                           "where there the mesh is allowed to move tangential to the "
                           "boundary. All tangential mesh movements along "
                           "those boundaries that have tangential material velocity "
                           "boundary conditions are allowed by default, this parameters "
                           "allows to generate mesh movements along other boundaries that are "
                           "open, or have prescribed material velocities or tractions."
                           "\n\n"
                           "The names of the boundaries listed here can either be "
                           "numbers (in which case they correspond to the numerical "
                           "boundary indicators assigned by the geometry object), or they "
                           "can correspond to any of the symbolic names the geometry object "
                           "may have provided for each part of the boundary. You may want "
                           "to compare this with the documentation of the geometry model you "
                           "use in your model.");
      }
      prm.leave_subsection ();

      std_cxx11::get<dim>(registered_plugins).declare_parameters (prm);
    }

    template <int dim>
    void FreeSurfaceHandler<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Free surface");
      {
        model_names
          = Utilities::split_string_list(prm.get("List of model names"));

        AssertThrow(Utilities::has_unique_entries(model_names),
                    ExcMessage("The list of strings for the parameter "
                               "'Initial temperature model/List of model names' contains entries more than once. "
                               "This is not allowed. Please check your parameter file."));

        // Create the list of tangential mesh movement boundary indicators
        try
          {
            const std::vector<types::boundary_id> x_additional_tangential_mesh_boundary_indicators
              = sim.geometry_model->translate_symbolic_boundary_names_to_ids(Utilities::split_string_list
                                                                             (prm.get ("Additional tangential mesh velocity boundary indicators")));

            tangential_mesh_boundary_indicators = sim.boundary_velocity_manager.get_tangential_boundary_velocity_indicators();
            tangential_mesh_boundary_indicators.insert(x_additional_tangential_mesh_boundary_indicators.begin(),
                                                       x_additional_tangential_mesh_boundary_indicators.end());
          }
        catch (const std::string &error)
          {
            AssertThrow (false, ExcMessage ("While parsing the entry <Free surface/Additional tangential "
                                            "mesh velocity boundary indicators>, there was an error. Specifically, "
                                            "the conversion function complained as follows: "
                                            + error));
          }
      }
      prm.leave_subsection ();

      // go through the list, create objects and let them parse
      // their own parameters
      for (unsigned int i=0; i<model_names.size(); ++i)
        {
          // create initial temperature objects
          mesh_deformation_objects.push_back (std_cxx11::shared_ptr<Interface<dim> >
                                              (std_cxx11::get<dim>(registered_plugins)
                                               .create_plugin (model_names[i],
                                                               "Initial temperature model::Model names")));

          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(&*mesh_deformation_objects.back()))
            sim->initialize_simulator (this->get_simulator());

          mesh_deformation_objects.back()->parse_parameters (prm);
          mesh_deformation_objects.back()->initialize ();
        }
    }



    template <int dim>
    void FreeSurfaceHandler<dim>::execute()
    {
      if (!sim.parameters.free_surface_enabled)
        return;

      TimerOutput::Scope timer (sim.computing_timer, "Mesh deformation");

      // Make the constraints for the elliptic problem.  On the free surface, we
      // constrain mesh velocity to be v.n, on free slip it is constrained to
      // be tangential, and on no slip boundaries it is zero.
      make_constraints();

      // Assemble and solve the vector Laplace problem which determines
      // the mesh displacements in the interior of the domain
      compute_mesh_displacements();

      // Interpolate the mesh velocity into the same
      // finite element space as used in the Stokes solve, which
      // is needed for the ALE corrections.
      interpolate_mesh_velocity();

      // After changing the mesh we need to rebuild things
      sim.rebuild_stokes_matrix = sim.rebuild_stokes_preconditioner = true;
    }



    template <int dim>
    void FreeSurfaceHandler<dim>::make_constraints()
    {
      if (!sim.parameters.free_surface_enabled)
        return;

      // Now construct the mesh displacement constraints
      mesh_displacement_constraints.clear();
      mesh_displacement_constraints.reinit(mesh_locally_relevant);

      // mesh_displacement_constraints can use the same hanging node
      // information that was used for mesh_vertex constraints.
      mesh_displacement_constraints.merge(mesh_vertex_constraints);

      // Add the vanilla periodic boundary constraints
      typedef std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int> > periodic_boundary_pairs;
      periodic_boundary_pairs pbp = sim.geometry_model->get_periodic_boundary_pairs();
      for (periodic_boundary_pairs::iterator p = pbp.begin(); p != pbp.end(); ++p)
        DoFTools::make_periodicity_constraints(free_surface_dof_handler, (*p).first.first, (*p).first.second, (*p).second, mesh_displacement_constraints);

      // Zero out the displacement for the zero-velocity boundary indicators
      for (std::set<types::boundary_id>::const_iterator p = sim.boundary_velocity_manager.get_zero_boundary_velocity_indicators().begin();
           p != sim.boundary_velocity_manager.get_zero_boundary_velocity_indicators().end(); ++p)
        VectorTools::interpolate_boundary_values (free_surface_dof_handler, *p,
                                                  ZeroFunction<dim>(dim), mesh_displacement_constraints);

      // Zero out the displacement for the prescribed velocity boundaries
      // if the boundary is not in the set of tangential mesh boundaries
      for (std::map<types::boundary_id, std::pair<std::string, std::vector<std::string> > >::const_iterator p = sim.boundary_velocity_manager.get_active_boundary_velocity_names().begin();
           p != sim.boundary_velocity_manager.get_active_boundary_velocity_names().end(); ++p)
        {
          if (tangential_mesh_boundary_indicators.find(p->first) == tangential_mesh_boundary_indicators.end())
            {
              VectorTools::interpolate_boundary_values (free_surface_dof_handler, p->first,
                                                        ZeroFunction<dim>(dim), mesh_displacement_constraints);
            }
        }

      sim.signals.pre_compute_no_normal_flux_constraints(sim.triangulation);
      // Make the no flux boundary constraints for boundaries with tangential mesh boundaries
      VectorTools::compute_no_normal_flux_constraints (free_surface_dof_handler,
                                                       /* first_vector_component= */
                                                       0,
                                                       tangential_mesh_boundary_indicators,
                                                       mesh_displacement_constraints, *sim.mapping);

      // make the periodic boundary indicators no displacement normal to the boundary
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
                                                       mesh_displacement_constraints, *sim.mapping);
      sim.signals.post_compute_no_normal_flux_constraints(sim.triangulation);

      // Ask all plugins to add their constraints
      // For the moment add constraints from all plugins into one matrix, then
      // merge that matrix with the existing constraints (respecting the existing
      // constraints as more important)
      ConstraintMatrix plugin_constraints(mesh_vertex_constraints.get_local_lines());

      for (unsigned int i=0; i<mesh_deformation_objects.size(); ++i)
        {
          ConstraintMatrix current_plugin_constraints(mesh_vertex_constraints.get_local_lines());

          mesh_deformation_objects[i]->deformation_constraints(free_surface_dof_handler,
                                                               current_plugin_constraints);

          const IndexSet local_lines = current_plugin_constraints.get_local_lines();
          for (auto index = local_lines.begin(); index != local_lines.end(); ++index)
            {
              if (current_plugin_constraints.is_constrained(*index))
                if (plugin_constraints.is_constrained(*index) == false)
                  {
                    plugin_constraints.add_line(*index);
                    plugin_constraints.set_inhomogeneity(*index, current_plugin_constraints.get_inhomogeneity(*index));
                  }
                else
                  {
                    const double inhomogeneity = plugin_constraints.get_inhomogeneity(*index);
                    plugin_constraints.set_inhomogeneity(*index, current_plugin_constraints.get_inhomogeneity(*index) + inhomogeneity);
                  }
            }
        }

      mesh_displacement_constraints.merge(plugin_constraints,ConstraintMatrix::left_object_wins);
      mesh_displacement_constraints.close();
    }



    template <int dim>
    void FreeSurfaceHandler<dim>::compute_mesh_displacements()
    {
      QGauss<dim> quadrature(free_surface_fe.degree + 1);
      UpdateFlags update_flags = UpdateFlags(update_values | update_JxW_values | update_gradients);
      FEValues<dim> fe_values (*sim.mapping, free_surface_fe, quadrature, update_flags);

      const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                         dofs_per_face = sim.finite_element.dofs_per_face,
                         n_q_points    = fe_values.n_quadrature_points;

      std::vector<types::global_dof_index> cell_dof_indices (dofs_per_cell);
      std::vector<unsigned int> face_dof_indices (dofs_per_face);
      Vector<double> cell_vector (dofs_per_cell);
      FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

      // We are just solving a Laplacian in each spatial direction, so
      // the degrees of freedom for different dimensions do not couple.
      Table<2,DoFTools::Coupling> coupling (dim, dim);
      coupling.fill(DoFTools::none);

      for (unsigned int c=0; c<dim; ++c)
        coupling[c][c] = DoFTools::always;

      LinearAlgebra::SparseMatrix mesh_matrix;
#ifdef ASPECT_USE_PETSC
      LinearAlgebra::DynamicSparsityPattern sp(mesh_locally_relevant);
#else
      TrilinosWrappers::SparsityPattern sp (mesh_locally_owned,
                                            mesh_locally_owned,
                                            mesh_locally_relevant,
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

      // carry out the solution
      FEValuesExtractors::Vector extract_vel(0);

      LinearAlgebra::Vector rhs, velocity_solution;
      rhs.reinit(mesh_locally_owned, sim.mpi_communicator);
      velocity_solution.reinit(mesh_locally_owned, sim.mpi_communicator);

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

      // Make the AMG preconditioner
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

      cg.solve (mesh_matrix, velocity_solution, rhs, preconditioner_stiffness);
      sim.pcout << "   Solving mesh velocity system... " << solver_control.last_step() <<" iterations."<< std::endl;

      mesh_displacement_constraints.distribute (velocity_solution);

      // Update the free surface mesh velocity vector
      fs_mesh_velocity = velocity_solution;

      // Update the mesh displacement vector
      LinearAlgebra::Vector distributed_mesh_displacements(mesh_locally_owned, sim.mpi_communicator);
      distributed_mesh_displacements = mesh_displacements;
      distributed_mesh_displacements.add(sim.time_step, velocity_solution);
      mesh_displacements = distributed_mesh_displacements;

    }


    template <int dim>
    void FreeSurfaceHandler<dim>::interpolate_mesh_velocity()
    {
      // Interpolate the mesh vertex velocity onto the Stokes velocity system for use in ALE corrections
      LinearAlgebra::BlockVector distributed_mesh_velocity;
      distributed_mesh_velocity.reinit(sim.introspection.index_sets.system_partitioning, sim.mpi_communicator);

      const std::vector<Point<dim> > support_points
        = sim.finite_element.base_element(sim.introspection.component_indices.velocities[0]).get_unit_support_points();

      Quadrature<dim> quad(support_points);
      UpdateFlags update_flags = UpdateFlags(update_values | update_JxW_values);
      FEValues<dim> fs_fe_values (*sim.mapping, free_surface_fe, quad, update_flags);
      FEValues<dim> fe_values (*sim.mapping, sim.finite_element, quad, update_flags);
      const unsigned int n_q_points = fe_values.n_quadrature_points,
                         dofs_per_cell = fe_values.dofs_per_cell;

      std::vector<types::global_dof_index> cell_dof_indices (dofs_per_cell);
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
            fs_fe_values[extract_vel].get_function_values(fs_mesh_velocity, velocity_values);
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
    void FreeSurfaceHandler<dim>::setup_dofs()
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

      mesh_displacements.reinit(mesh_locally_owned, mesh_locally_relevant, sim.mpi_communicator);
      fs_mesh_velocity.reinit(mesh_locally_owned, mesh_locally_relevant, sim.mpi_communicator);

      // if we are just starting, we need to initialize the mesh displacement vector.
      if (sim.timestep_number == 0)
        mesh_displacements = 0.;

      // We would like to make sure that the mesh stays conforming upon
      // redistribution, so we construct mesh_vertex_constraints, which
      // keeps track of hanging node constraints.
      // Note: this would be a more natural fit in make_constraints(),
      // but we would like to be able to apply vertex constraints directly
      // after setup_dofs(), as is done, for instance, during mesh
      // refinement.
      mesh_vertex_constraints.clear();
      mesh_vertex_constraints.reinit(mesh_locally_relevant);

      DoFTools::make_hanging_node_constraints(free_surface_dof_handler, mesh_vertex_constraints);

      // We can safely close this now
      mesh_vertex_constraints.close();

      // Now reset the mapping of the simulator to be something that captures mesh deformation in time.
      sim.mapping.reset (new MappingQ1Eulerian<dim, LinearAlgebra::Vector> (free_surface_dof_handler,
                                                                            mesh_displacements));
    }



    template <int dim>
    const std::vector<std::string> &
    FreeSurfaceHandler<dim>::get_active_mesh_deformation_names () const
    {
      return model_names;
    }



    template <int dim>
    const std::vector<std_cxx11::shared_ptr<Interface<dim> > > &
    FreeSurfaceHandler<dim>::get_active_mesh_deformation_models () const
    {
      return mesh_deformation_objects;
    }



    template <int dim>
    std::string
    get_valid_model_names_pattern ()
    {
      return std_cxx11::get<dim>(registered_plugins).get_pattern_of_names ();
    }



    template <int dim>
    void
    FreeSurfaceHandler<dim>::write_plugin_graph (std::ostream &out)
    {
      std_cxx11::get<dim>(registered_plugins).write_plugin_graph ("Mesh deformation interface",
                                                                  out);
    }
  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace internal
  {
    namespace Plugins
    {
      template <>
      std::list<internal::Plugins::PluginList<MeshDeformation::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<MeshDeformation::Interface<2> >::plugins = 0;
      template <>
      std::list<internal::Plugins::PluginList<MeshDeformation::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<MeshDeformation::Interface<3> >::plugins = 0;
    }
  }

  namespace MeshDeformation
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  template class FreeSurfaceHandler<dim>; \
  \
  template \
  std::string \
  get_valid_model_names_pattern<dim> ();

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}
