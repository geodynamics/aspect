/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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


#include <iostream>

#include <aspect/global.h>
#include <aspect/mesh_deformation/fastscapecc.h>
#include <aspect/mesh_deformation/fastscapecc_adapter.h>
#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/vector_tools_interpolate.h>
#include <deal.II/lac/affine_constraints.h>
#include <aspect/simulator.h>
#include <aspect/simulator/assemblers/interface.h>
#include <aspect/simulator_signals.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/fe_field_function.h>
#include <typeinfo>

namespace aspect
{
  namespace MeshDeformation
  {
    template <int dim>
    FastScapecc<dim>::FastScapecc()
      // : surface_mesh_dof_handler(surface_mesh) // Link DoFHandler to surface_mesh
    {}

    template <int dim>
    void FastScapecc<dim>::initialize()
    {
      const auto &geom_model = this->get_geometry_model();

      this->get_pcout() << "Geometry model type: " << typeid(geom_model).name() << std::endl;

      // Runtime dispatch to correct overload
      if (const auto *box = dynamic_cast<const GeometryModel::Box<dim> *>(&geom_model))
        init_surface_mesh(*box);
      else if (const auto *sphere = dynamic_cast<const GeometryModel::SphericalShell<dim> *>(&geom_model))
        init_surface_mesh(*sphere);
      else
        AssertThrow(false, ExcMessage("FastScapecc plugin only supports Box or Spherical Shell geometries."));

      n_grid_nodes = surface_mesh.n_active_cells();
    }

    template <int dim>
    void FastScapecc<dim>::init_surface_mesh(const GeometryModel::Box<dim> &geom_model)
    {
      this->get_pcout() << "Box geometry detected. Initializing FastScape for Box geometry..." << std::endl;

      const auto origin = geom_model.get_origin();
      const auto extent = geom_model.get_extents();

      const Point<dim - 1> p1(origin[0], origin[1]);
      const Point<dim - 1> p2(origin[0] + extent[0], origin[1] + extent[1]);

      auto rep = geom_model.get_repetitions();
      std::vector<unsigned int> repetitions(rep.begin(), rep.begin() + (dim - 1));

      // Store for indexing
      grid_extent[0] = std::make_pair(origin[0], origin[0] + extent[0]);
      grid_extent[1] = std::make_pair(origin[1], origin[1] + extent[1]);
      nx = repetitions[0];
      ny = repetitions[1];
      dx = (grid_extent[0].second - grid_extent[0].first) / static_cast<double>(nx);
      dy = (grid_extent[1].second - grid_extent[1].first) / static_cast<double>(ny);

      // Generate surface mesh
      GridGenerator::subdivided_hyper_rectangle(surface_mesh, repetitions, p1, p2);

      surface_fe = std::make_shared<FE_Q<dim - 1, dim>>(1);
      surface_mesh_dof_handler.reinit(surface_mesh);
      surface_mesh_dof_handler.distribute_dofs(*surface_fe);

      surface_constraints.clear();
      DoFTools::make_hanging_node_constraints(surface_mesh_dof_handler, surface_constraints);
      surface_constraints.close();
}



    template <int dim>
    void FastScapecc<dim>::init_surface_mesh(const GeometryModel::SphericalShell<dim> &geom_model)
    {
      this->get_pcout() << "Spherical Shell geometry detected. Initializing FastScape for Spherical Shell geometry..." << std::endl;

      const Point<dim> center;
      GridGenerator::hyper_sphere(surface_mesh, center, geom_model.outer_radius());
      surface_mesh.refine_global(3);


      surface_fe = std::make_shared<FE_Q<dim - 1, dim>>(1);
      surface_mesh_dof_handler.reinit(surface_mesh);

      surface_mesh_dof_handler.distribute_dofs(*surface_fe);

      surface_constraints.clear();
      DoFTools::make_hanging_node_constraints(surface_mesh_dof_handler, surface_constraints);
      surface_constraints.close();
    }


    template <int dim>
    void FastScapecc<dim>::project_surface_solution(const std::set<types::boundary_id> &boundary_ids)
    {
      TimerOutput::Scope timer_section(this->get_computing_timer(), "Project surface solution");

      // Define the quadrature rule for projection
      QGauss<2> quadrature(surface_mesh_dof_handler.get_fe().degree + 1);

      // Initialize the surface solution vector
      IndexSet locally_relevant_dofs_surface;
      DoFTools::extract_locally_relevant_dofs(surface_mesh_dof_handler, locally_relevant_dofs_surface);
      surface_solution.reinit(surface_mesh_dof_handler.locally_owned_dofs(),
                              locally_relevant_dofs_surface,
                              this->get_mpi_communicator());

      // Initialize the boundary solution vector
      IndexSet locally_relevant_dofs_main;
      DoFTools::extract_locally_relevant_dofs(this->get_dof_handler(), locally_relevant_dofs_main);
      boundary_solution.reinit(this->get_dof_handler().locally_owned_dofs(),
                               locally_relevant_dofs_main,
                               this->get_mpi_communicator());

      const types::boundary_id relevant_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
                {
                  if (cell->face(face)->at_boundary() &&
                      cell->face(face)->boundary_id() == relevant_boundary)
                    {
                      for (unsigned int vertex = 0; vertex < GeometryInfo<dim - 1>::vertices_per_face; ++vertex)
                        {
                          const unsigned int dof_index = cell->face(face)->vertex_dof_index(vertex, 0);
                          boundary_solution[dof_index] = this->get_solution()[dof_index];
                        }
                    }
                }
            }
        }

      boundary_solution.compress(VectorOperation::insert);

    }


    template <int dim>
    unsigned int FastScapecc<dim>::vertex_index(const Point<dim> &p) const
    {
      if constexpr (dim == 3)
        {
          const unsigned int i = static_cast<unsigned int>((p[0] - grid_extent[0].first) / dx + 0.5);
          const unsigned int j = static_cast<unsigned int>((p[1] - grid_extent[1].first) / dy + 0.5);
          return j * nx + i;
        }
      else if constexpr (dim == 2)
        {
          const unsigned int i = static_cast<unsigned int>((p[0] - grid_extent[0].first) / dx + 0.5);
          return i;
        }
      else
        AssertThrow(false, ExcMessage("Unsupported dimension in vertex_index()."));
    }


    template <int dim>
    void FastScapecc<dim>::compute_velocity_constraints_on_boundary(
        const DoFHandler<dim> &mesh_deformation_dof_handler,
        AffineConstraints<double> &mesh_velocity_constraints,
        const std::set<types::boundary_id> &boundary_ids) const
    {
      if (this->get_timestep_number() == 0)
        return;

      TimerOutput::Scope timer_section(this->get_computing_timer(), "FastScape plugin");

      const_cast<FastScapecc<dim> *>(this)->project_surface_solution(boundary_ids);

      // Identify geometry model type
      const auto *spherical_model = dynamic_cast<const GeometryModel::SphericalShell<dim> *>(&this->get_geometry_model());
      const auto *box_model = dynamic_cast<const GeometryModel::Box<dim> *>(&this->get_geometry_model());

      AssertThrow(spherical_model || box_model, ExcMessage("Unsupported geometry model in FastScapecc."));

      // Temporary storage: elevation, vertex_index, and vertical/radial velocity
      std::vector<std::vector<double>> temporary_variables(dim + 2, std::vector<double>());

      for (const auto &cell : surface_mesh_dof_handler.active_cell_iterators())
        {
          for (unsigned int v = 0; v < GeometryInfo<dim - 1>::vertices_per_cell; ++v)
            {
              const Point<dim> vertex = cell->vertex(v);
              const unsigned int dof_index = cell->vertex_dof_index(v, 0);
              const double surface_value = surface_solution[dof_index];

              // Compute elevation
              double elevation = vertex[2];
              if (spherical_model)
                elevation -= spherical_model->outer_radius();

              // Compute velocity (radial or vertical)
              double velocity = 0.0;
              if (spherical_model)
                velocity = (vertex[0] * surface_value +
                            vertex[1] * surface_value +
                            vertex[2] * surface_value) / vertex.norm();
              else if (box_model)
                velocity = surface_value;

              // Get 1D or 2D grid index from physical position
              const unsigned int index = this->vertex_index(vertex);

              // Fill temporary variables
              temporary_variables[0].push_back(elevation);
              temporary_variables[1].push_back(static_cast<double>(index));
              temporary_variables[2].push_back(velocity * year_in_seconds);
            }
        }
    
      std::vector<double> V(n_grid_nodes);

      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          // Initialize the variables that will be sent to FastScape.
          std::vector<double> h(n_grid_nodes, std::numeric_limits<double>::max());
          std::vector<double> vz(n_grid_nodes);
          std::vector<double> vy(n_grid_nodes);
          std::vector<double> h_old(n_grid_nodes);

          for (unsigned int i = 0; i < temporary_variables[1].size(); ++i)
            {
              int index = static_cast<int>(temporary_variables[1][i]);
              h[index] = temporary_variables[0][i];
              vz[index] = temporary_variables[2][i];
              vy[index] = temporary_variables[3][i];
            }

          for (unsigned int p = 1; p < Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++p)
            {
              MPI_Status status;
              MPI_Probe(p, 42, this->get_mpi_communicator(), &status);
              int incoming_size = 0;
              MPI_Get_count(&status, MPI_DOUBLE, &incoming_size);

              for (unsigned int i = 0; i < temporary_variables.size(); ++i)
                {
                  temporary_variables[i].resize(incoming_size);
                }

              for (unsigned int i = 0; i < temporary_variables.size(); ++i)
                MPI_Recv(&temporary_variables[i][0], incoming_size, MPI_DOUBLE, p, 42, this->get_mpi_communicator(), &status);

              for (unsigned int i = 0; i < temporary_variables[1].size(); ++i)
                {
                  int index = static_cast<int>(temporary_variables[1][i]);
                  h[index] = temporary_variables[0][i];
                  vz[index] = temporary_variables[dim+1][i];
                }
            }

          for (unsigned int i = 0; i < n_grid_nodes; ++i)
            {
              h_old[i] = h[i];
            }

          auto elevation = xt::adapt(h);
          auto elevation_old = xt::adapt(h_old);

          const double aspect_timestep_in_years = this->get_timestep() / year_in_seconds;

          unsigned int fastscape_iterations = fastscape_steps_per_aspect_step;
          double fastscape_timestep_in_years = aspect_timestep_in_years / fastscape_iterations;
          while (fastscape_timestep_in_years > maximum_fastscape_timestep)
            {
              fastscape_iterations *= 2;
              fastscape_timestep_in_years *= 0.5;
            }

          double cell_area = surface_mesh.begin_active()->measure();

          // 1. Create the FastScape grid adapter
          auto grid = GridAdapterType(const_cast<SurfaceMeshType &>(surface_mesh), cell_area);

          // 2. Define the node status array BEFORE creating the flow graph //Benoit
          auto node_status_array = xt::zeros<fastscapelib::node_status>({n_grid_nodes});
          grid.set_nodes_status(node_status_array);

          // 3. Create the flow graph 
          auto flow_graph = FlowGraphType(
              grid,
              { fastscapelib::single_flow_router() }
          );

          // Pick an arbitrary base node for the flow to be computed
          std::vector<std::size_t> base_level_nodes = {0};
          // Set base level nodes
          flow_graph.set_base_levels(base_level_nodes);
          // 4. Now you can build the erosion model
          auto spl_eroder = fastscapelib::spl_eroder<FlowGraphType>(
              fastscapelib::make_spl_eroder(flow_graph, kff, n, m, kdd)
          );


          // Create data arrays using grid->shape()
          auto uplift_rate = xt::adapt(vz);
          xt::xarray<double> drainage_area = xt::zeros<double>(grid.shape());
          xt::xarray<double> sediment_flux = xt::zeros<double>(grid.shape());

          // Start erosion loop
          xt::xarray<double> uplifted_elevation = elevation + fastscape_timestep_in_years * uplift_rate;

          for (unsigned int i = 0; i < fastscape_iterations; ++i)
          {
            uplifted_elevation = elevation + fastscape_timestep_in_years * uplift_rate;
            flow_graph.update_routes(uplifted_elevation);
            flow_graph.accumulate(drainage_area, 1.0);
            auto spl_erosion = spl_eroder.erode(uplifted_elevation, drainage_area, fastscape_timestep_in_years);
            sediment_flux = flow_graph.accumulate(spl_erosion);
            elevation = uplifted_elevation - spl_erosion;
          }

          // Compute erosion velocities
          for (unsigned int i = 0; i < n_grid_nodes; ++i)
            V[i] = (elevation[i] - elevation_old[i]) / aspect_timestep_in_years;

          // Broadcast V to all processes
          MPI_Bcast(&V[0], n_grid_nodes, MPI_DOUBLE, 0, this->get_mpi_communicator());
        }
      else
        {
          for (unsigned int i = 0; i < temporary_variables.size(); ++i)
            MPI_Ssend(&temporary_variables[i][0], temporary_variables[1].size(), MPI_DOUBLE, 0, 42, this->get_mpi_communicator());
          MPI_Bcast(&V[0], n_grid_nodes, MPI_DOUBLE, 0, this->get_mpi_communicator());
        }

      // Step 1: Interpolation function from position to velocity
      auto erosion_function = [&](const Point<dim> &p) -> double
      {
      // Closest-point search (replace with a spatial tree for performance)
      double min_dist = std::numeric_limits<double>::max();
      unsigned int best_index = 0;

      unsigned int i = 0;
      for (const auto &cell : surface_mesh.active_cell_iterators())
        for (unsigned int v = 0; v < GeometryInfo<dim-1>::vertices_per_cell; ++v, ++i)
          {
            const Point<dim> &node = cell->vertex(v);
            double dist = node.distance(p);
            if (dist < min_dist)
              {
                min_dist = dist;
                best_index = i;
              }
          }

      return V[best_index];
      };

      VectorFunctionFromScalarFunctionObject<dim> radial_velocity_field(
          erosion_function,
          dim - 1, // project onto radial component
          dim      // full space dimension
      );


      for (const auto boundary_id : boundary_ids)
        VectorTools::interpolate_boundary_values(
            this->get_mapping(),
            mesh_deformation_dof_handler,
            boundary_id,
            radial_velocity_field,
            mesh_velocity_constraints);

      this->get_pcout() << "Applying erosion velocities to mesh boundary." << std::endl;

    }

    template <int dim>
    bool
    FastScapecc<dim>::
    needs_surface_stabilization () const
    {
      return true;
    }


    template <int dim>
    void FastScapecc<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection ("Fastscapecc");
        {
          prm.declare_entry("Number of steps", "10",
                            Patterns::Integer(),
                            "Number of steps per ASPECT timestep");
          prm.declare_entry("Maximum timestep", "10e3",
                            Patterns::Double(0),
                            "Maximum timestep for FastScape. Units: $\\{yrs}$");
          prm.declare_entry("Additional fastscape refinement levels", "0",
                            Patterns::Integer(),
                            "How many levels above ASPECT FastScape should be refined.");
          prm.declare_entry("Fastscape seed", "1000",
                            Patterns::Integer(),
                            "Seed used for adding an initial noise to FastScape topography based on the initial noise magnitude.");
          prm.declare_entry("Maximum surface refinement level", "1",
                            Patterns::Integer(),
                            "This should be set to the highest ASPECT refinement level expected at the surface.");
          prm.declare_entry("Surface refinement difference", "0",
                            Patterns::Integer(),
                            "The difference between the lowest and highest refinement level at the surface. E.g., if three resolution "
                            "levels are expected, this would be set to 2.");
          prm.declare_entry ("Use velocities", "true",
                             Patterns::Bool (),
                             "Flag to use FastScape advection and uplift.");
          prm.declare_entry("Precision", "0.001",
                            Patterns::Double(),
                            "Precision value for how close a ASPECT node must be to the FastScape node for the value to be transferred.");
          prm.declare_entry("Initial noise magnitude", "5",
                            Patterns::Double(),
                            "Maximum topography change from the initial noise. Units: $\\{m}$");

          prm.enter_subsection ("Erosional parameters");
          {
            prm.declare_entry("Drainage area exponent", "0.4",
                              Patterns::Double(),
                              "Exponent for drainage area.");
            prm.declare_entry("Slope exponent", "1",
                              Patterns::Double(),
                              "The  slope  exponent  for  SPL (n).  Generally  m/n  should  equal  approximately 0.4");
            prm.declare_entry("Bedrock river incision rate", "1e-5",
                              Patterns::Double(),
                              "River incision rate for bedrock in the Stream Power Law. Units: $\\{m^(1-2*drainage_area_exponent)/yr}$");
            prm.declare_entry("Sediment river incision rate", "-1",
                              Patterns::Double(),
                              "River incision rate for sediment in the Stream Power Law. -1 sets this to the bedrock river incision rate. Units: $\\{m^(1-2*drainage_area_exponent)/yr}$ ");
            prm.declare_entry("Bedrock diffusivity", "1e-2",
                              Patterns::Double(),
                              "Transport coefficient (diffusivity) for bedrock. Units: $\\{m^2/yr}$ ");
            prm.declare_entry("Sediment diffusivity", "-1",
                              Patterns::Double(),
                              "Transport coefficient (diffusivity) for sediment. -1 sets this to the bedrock diffusivity. Units: $\\{m^2/yr}$");
            prm.declare_entry("Elevation factor", "1",
                              Patterns::Double(),
                              "Amount to multiply kf and kd by past given orographic elevation control.");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void FastScapecc<dim>::parse_parameters(ParameterHandler &prm)
    {
      end_time = prm.get_double ("End time");
      if (prm.get_bool ("Use years in output instead of seconds") == true)
        end_time *= year_in_seconds;

      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection("Fastscapecc");
        {
          fastscape_steps_per_aspect_step = prm.get_integer("Number of steps");
          maximum_fastscape_timestep = prm.get_double("Maximum timestep");
          additional_refinement_levels = prm.get_integer("Additional fastscape refinement levels");
          fs_seed = prm.get_integer("Fastscape seed");
          maximum_surface_refinement_level = prm.get_integer("Maximum surface refinement level");
          surface_refinement_difference = prm.get_integer("Surface refinement difference");
          precision = prm.get_double("Precision");
          noise_h = prm.get_double("Initial noise magnitude");

          if (!this->convert_output_to_years())
            {
              maximum_fastscape_timestep /= year_in_seconds;
            }

          prm.enter_subsection("Erosional parameters");
          {
            m = prm.get_double("Drainage area exponent");
            n = prm.get_double("Slope exponent");
            kfsed = prm.get_double("Sediment river incision rate");
            kff = prm.get_double("Bedrock river incision rate");
            kdsed = prm.get_double("Sediment diffusivity");
            kdd = prm.get_double("Bedrock diffusivity");

            // if (!this->convert_output_to_years())
            //   {
            //     kff *= year_in_seconds;
            //     kdd *= year_in_seconds;
            //     kfsed *= year_in_seconds;
            //     kdsed *= year_in_seconds;
            //   }

          }
          prm.leave_subsection();

        }
        prm.leave_subsection();
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
    ASPECT_REGISTER_MESH_DEFORMATION_MODEL(FastScapecc,
                                           "fastscapecc",
                                           "A plugin, which prescribes the surface mesh to "
                                           "deform according to an analytically prescribed "
                                           "function. Note that the function prescribes a "
                                           "deformation velocity, i.e. the return value of "
                                           "this plugin is later multiplied by the time step length "
                                           "to compute the displacement increment in this time step. "
                                           "The format of the "
                                           "functions follows the syntax understood by the "
                                           "muparser library, see Section~\\ref{sec:muparser-format}.")
  }
}
