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

#ifdef ASPECT_WITH_FASTSCAPELIB

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
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/base/mpi.h>  // For Utilities::MPI::this_mpi_process
#include <aspect/utilities.h>  // For Utilities::create_directory
#include <filesystem>          // C++17 or later


namespace aspect
{
  namespace MeshDeformation
  {
    template <int dim>
    FastScapecc<dim>::FastScapecc()
      :
      surface_mesh(MPI_COMM_WORLD),
      surface_fe(1)  // Use a Q1 element for the surface mesh
      // , surface_mesh_dof_handler(surface_mesh) // Link DoFHandler to surface_mesh
    {}

    template <int dim>
    void FastScapecc<dim>::initialize()
    {
      const auto &geom_model = this->get_geometry_model();
      init_surface_mesh(geom_model);

      //Fastscapelib operates on grid nodes, each of which conceptually represents a cell center
      //and has associated area and neighbor connections.
      // n_grid_nodes = surface_mesh.n_active_cells();
      n_grid_nodes = surface_mesh.n_used_vertices();
    }

    template <int dim>
    void FastScapecc<dim>::init_surface_mesh(const GeometryModel::Interface<dim> &geom_model)
    {
      // First create the surface triangulation objects. This is necessarily different
      // between geometries.
      if (const auto *box = dynamic_cast<const GeometryModel::Box<dim> *>(&geom_model))
        {
          this->get_pcout() << "Box geometry detected. Initializing FastScape for Box geometry..." << std::endl;

          // Get origin and extent of the ASPECT domain
          const auto origin = box->get_origin();
          const auto extent = box->get_extents();

          for (unsigned int i = 0; i < dim; ++i)
            {
              box_grid_extent[i].first = origin[i];
              box_grid_extent[i].second = origin[i] + extent[i];
            }

          // Build surface mesh corners from grid extent
          Point<dim - 1> p1, p2;
          for (unsigned int i = 0; i < dim - 1; ++i)
            {
              p1[i] = box_grid_extent[i].first;
              p2[i] = box_grid_extent[i].second;
            }

          // Extract surface repetitions from full domain repetitions
          std::vector<unsigned int> surface_repetitions(repetitions.begin(), repetitions.begin() + (dim - 1));

          // Apply surface refinement factor
          const unsigned int total_refinement = surface_refinement_level + additional_refinement_levels;
          const unsigned int refinement_factor = static_cast<unsigned int>(std::pow(2.0, total_refinement));

          for (unsigned int i = 0; i < dim - 1; ++i)
            surface_repetitions[i] *= refinement_factor;

          // Store grid dimensions for indexing and spacing
          nx = surface_repetitions[0]+1;
          ny = surface_repetitions[1]+1;
          dx = (box_grid_extent[0].second - box_grid_extent[0].first) / static_cast<double>(nx-1);
          dy = (box_grid_extent[1].second - box_grid_extent[1].first) / static_cast<double>(ny-1);

          // Create refined surface mesh
          GridGenerator::subdivided_hyper_rectangle(surface_mesh,
                                                    surface_repetitions,
                                                    p1,
                                                    p2,
                                                    /* colorize */ true);
          Assert (surface_mesh.n_used_vertices() == nx*ny, ExcInternalError());
        }
      else if (const auto *spherical_shell = dynamic_cast<const GeometryModel::SphericalShell<dim> *>(&geom_model))
        {
          this->get_pcout() << "Spherical Shell geometry detected. Initializing FastScape for Spherical Shell geometry..." << std::endl;

          const Point<dim> center;
          GridGenerator::hyper_sphere(surface_mesh, center, spherical_shell->outer_radius());
          surface_mesh.refine_global(3); // Adjust as needed

          // Create a unique index for each surface vertex
          unsigned int counter = 0;
          for (const auto &cell : surface_mesh.active_cell_iterators())
            for (unsigned int v = 0; v < GeometryInfo<dim - 1>::vertices_per_cell; ++v)
              {
                const Point<dim> &vertex = cell->vertex(v);
                if (spherical_vertex_index_map.find(vertex) == spherical_vertex_index_map.end())
                  {
                    spherical_vertex_index_map[vertex] = counter++;
                  }
              }

          n_grid_nodes = spherical_vertex_index_map.size();
        }
      else
        AssertThrow(false, ExcMessage("FastScapecc plugin only supports Box or Spherical Shell geometries."));

      // Having create the mesh, now set up the DoFHandlers and constraints objects
      surface_mesh_dof_handler.reinit(surface_mesh);
      surface_mesh_dof_handler.distribute_dofs(surface_fe);

      surface_constraints.clear();
      DoFTools::make_hanging_node_constraints(surface_mesh_dof_handler, surface_constraints);
      surface_constraints.close();
    }



    template <int dim>
    void FastScapecc<dim>::project_surface_solution(const std::set<types::boundary_id> & /*boundary_ids*/,
                                                    dealii::LinearAlgebra::distributed::Vector<double> &surface_vertical_velocity,
                                                    dealii::LinearAlgebra::distributed::Vector<double> &surface_elevation ) const
    {
      TimerOutput::Scope timer_section(this->get_computing_timer(), "Project surface solution");

      const auto *spherical_model = dynamic_cast<const GeometryModel::SphericalShell<dim> *>(&this->get_geometry_model());
      const auto *box_model       = dynamic_cast<const GeometryModel::Box<dim> *>(&this->get_geometry_model());

      AssertThrow(spherical_model || box_model,
                  ExcMessage("FastScapecc only supports Box or SphericalShell geometries."));

      // Initialize the surface solution vector on the FastScape surface mesh
      const IndexSet locally_relevant_dofs_surface
        = DoFTools::extract_locally_relevant_dofs(surface_mesh_dof_handler);

      surface_vertical_velocity.reinit(surface_mesh_dof_handler.locally_owned_dofs(),
                                       locally_relevant_dofs_surface,
                                       this->get_mpi_communicator());

      surface_elevation.reinit(surface_mesh_dof_handler.locally_owned_dofs(),
                               locally_relevant_dofs_surface,
                               this->get_mpi_communicator());
      const types::boundary_id top_boundary =
        this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");
      // std::cout<<"here it works 0a"<<std::endl;

      // Use `volume_cell` for volume mesh iterations,
      // and `surface_cell` for surface mesh iterations,

      for (const auto &volume_cell : this->get_dof_handler().active_cell_iterators())
        if (volume_cell->is_locally_owned())
          {
            for (unsigned int face = 0; face < volume_cell->n_faces(); ++face)
              {
                if (!volume_cell->face(face)->at_boundary() ||
                    volume_cell->face(face)->boundary_id() != top_boundary)
                  continue;

                for (unsigned int v = 0; v < volume_cell->face(face)->n_vertices(); ++v)
                  {
                    // std::cout<<"here it works 0a"<<std::endl;

                    const Point<dim> pos = volume_cell->face(face)->vertex(v);
                    // std::cout<<"here it works 0b"<<std::endl;

                    // Extract the full velocity vector at the vertex
                    Tensor<1, dim> velocity;
                    for (unsigned int d = 0; d < dim; ++d)
                      {
                        // Query velocity Dofs using introspection
                        const unsigned int component_index = this->introspection().component_indices.velocities[d];
                        velocity[d] = this->get_solution()[volume_cell->face(face)->vertex_dof_index(v, component_index)];
                      }
                    // std::cout<<"here it works 0c"<<std::endl;


                    // Use the gravity model to obtain a "vertical" direction
                    const Tensor<1, dim> gravity = this->get_gravity_model().gravity_vector(pos);
                    Tensor<1, dim> gravity_dir = gravity;

                    const double gravity_norm = gravity.norm();
                    if (gravity_norm > 0.0)
                      gravity_dir /= gravity_norm;
                    // std::cout<<"here it works 0d"<<std::endl;


                    const double vertical_velocity = velocity * gravity_dir;
                    const double height = this->get_geometry_model().height_above_reference_surface(pos);
                    const unsigned int index = this->vertex_index(pos);
                    // std::cout<<"here it works 0e"<<std::endl;
                    // std::cout<< "Surface pos: " << pos <<std::endl;
                    // std::cout   << ", Index: " << index << std::endl;
                    // std::cout<< "Vertical velocity: " << vertical_velocity<< std::endl;
                    // std::cout  << ", velocity: " << surface_vertical_velocity[index]<<std::endl;
                    // std::cout   << ", topography: " << surface_elevation[index] <<std::endl;

                    surface_vertical_velocity[index] = vertical_velocity;
                    // std::cout<<"here it works 0f"<<std::endl;

                    surface_elevation[index] = height;
                    // std::cout<<"here it works 0g"<<std::endl;

                  }
              }
          }
      surface_vertical_velocity.compress(dealii::VectorOperation::insert);
      surface_elevation.compress(dealii::VectorOperation::insert);

      this->get_pcout() << "Projected surface solution onto FastScape mesh." << std::endl;

    }

    template <int dim>
    unsigned int FastScapecc<dim>::vertex_index(const Point<dim> &p) const
    {
      if (const auto *box_model = dynamic_cast<const GeometryModel::Box<dim> *>(&this->get_geometry_model()))
        {
          if constexpr (dim == 3)
            {
              const unsigned int i =std::min(static_cast<unsigned int>((p[0] - box_grid_extent[0].first) / dx + 0.5), nx - 1);
              const unsigned int j =std::min(static_cast<unsigned int>((p[1] - box_grid_extent[1].first) / dy + 0.5), ny - 1);
              const unsigned int index = j * nx + i;

              // std::cout << "[DEBUG vertex_index] p: (" << x << ", " << y << "), "
              //           << "box_grid_extent x0: " << x0 << ", y0: " << y0 << ", "
              //           << "dx: " << dx << ", dy: " << dy << ", "
              //           << "i: " << i << ", j: " << j << ", nx: " << nx << ", "
              //           << "index: " << index << std::endl;

              return index;
            }
          else if constexpr (dim == 2)
            {
              const double x = p[0];
              const double x0 = box_grid_extent[0].first;
              const unsigned int i = static_cast<unsigned int>((x - x0) / dx + 0.5);

              std::cout << "[DEBUG vertex_index 2D] p: " << x
                        << ", box_grid_extent x0: " << x0
                        << ", dx: " << dx << ", i: " << i << std::endl;

              return i;
            }
        }
      else if (const auto *spherical_model = dynamic_cast<const GeometryModel::SphericalShell<dim> *>(&this->get_geometry_model()))
        {
          const auto it = spherical_vertex_index_map.find(p);
          AssertThrow(it != spherical_vertex_index_map.end(),
                      ExcMessage("Point not found in spherical_vertex_index_map."));
          return it->second;
        }
      else
        {
          AssertThrow(false, ExcMessage("Unsupported geometry in vertex_index()."));
          return numbers::invalid_unsigned_int;
        }
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

      // Step 1: Project current surface velocity from ASPECT solution
      dealii::LinearAlgebra::distributed::Vector<double> surface_vertical_velocity, surface_elevation;
      // std::cout<<"here it works 0"<<std::endl;
      this->project_surface_solution(boundary_ids, surface_vertical_velocity, surface_elevation);
      // std::cout<<"here it works 1"<<std::endl;

      const auto *spherical_model = dynamic_cast<const GeometryModel::SphericalShell<dim> *>(&this->get_geometry_model());
      const auto *box_model = dynamic_cast<const GeometryModel::Box<dim> *>(&this->get_geometry_model());
      AssertThrow(spherical_model || box_model, ExcMessage("Unsupported geometry model in FastScapecc."));

// TODO: You are using these three arrays in the following way:
//   temp[0][i] = topography value at some location
//   temp[1][i] = DoF index within surface_dof_handler for this location
//   temp[2][i] = surface_uplift_rate for this location
// I think this would be more transparent if you had the following:
// struct ValuesAtSurfaceVertex
// {
//   double topography;
//   double surface_uplift_rate;
// };
// and then you store everything in a
// std::multimap<types::global_dof_index, ValuesAtSurfaceVertex> surface_vertex_to_surface_values;
// Note that here I chose a multimap because you visit vertices more than once from all of
// the adjacent cells, at least in the current code you enter the same vertex into
//the temp vectors multiple times. I *think* that you should get the same values every
//time you visit a vertex, so perhaps a std::map is a better choice.
      std::vector<std::vector<double>> temporary_variables(3, std::vector<double>());

// TODO: Based on the convention, one should then rename cell -> surface_cell
      for (const auto &surface_cell : surface_mesh_dof_handler.active_cell_iterators())
        {
          for (unsigned int vertex_index = 0; vertex_index < surface_cell->n_vertices(); ++vertex_index)
            {
              const Point<dim> vertex = surface_cell->vertex(vertex_index);
              const unsigned int dof_index = surface_cell->vertex_dof_index(vertex_index, 0);

              const double surface_uplift_rate = surface_vertical_velocity[dof_index];
              const double topography = surface_elevation[dof_index];

              const unsigned int index = this->vertex_index(vertex);

              // Optional debug output
              if (this->get_timestep_number() == 1)
                this->get_pcout() << "Vertex at z=" << vertex[dim - 1]
                                  << ", DoF index: " << dof_index
                                  << ", vertex_index: " << index
                                  << ", uplift rate: " << surface_uplift_rate << " m/s"
                                  << ", topography: " << topography << " m"
                                  << std::endl;

              temporary_variables[0].push_back(topography);
              temporary_variables[1].push_back(static_cast<double>(index));
              temporary_variables[2].push_back(surface_uplift_rate * year_in_seconds);  // convert to m/year
            }
        }


      std::vector<double> V(n_grid_nodes);

// TODO: If you do the work in the following only of process 0, don't the other processes
// have to send their data to process 0? I.e., don't we need a 'gather' operation here?
//
// Update: On reading further, I see that you actually do that. Let's move the 'else' branch
// here, doing something like
//    if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) != 0)
//      send data to process zero
// This makes it easier to see that in order for this to work, we first have to send the
// data before we receive it in the ==0 block.
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          // Initialize the variables that will be sent to FastScape.
          std::vector<double> h(n_grid_nodes, std::numeric_limits<double>::max());
          std::vector<double> vz(n_grid_nodes);

// TODO: This will then simply be a loop over all elements of the map/multimap above.
// Note that right now you encounter the same vertex multiple times, and you simply
// overwrite early content with later content (which is ok).
          // First copy our own data into the h and vz arrays:
          for (unsigned int i = 0; i < temporary_variables[1].size(); ++i)
            {
              int index = static_cast<int>(temporary_variables[1][i]);
              h[index] = temporary_variables[0][i];
              vz[index] = temporary_variables[2][i];
            }

          // Then also retrieve what the other processes sent to us and put
          // them into the same arrays as well:
          for (unsigned int p = 1; p < Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++p)
            {
// TODO: In order to send/receive the std::map or std::multimap mentioned above, you can use
// Utilities::MPI::isend/irecv. These are easier to use than what you have here because
// they work on *any* data structure.

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
                  vz[index] = temporary_variables[2][i];
                }
            }

          // Save the current elevation so that we can later take the difference to
          // the one after calling FastScape:
          const std::vector<double> h_old = h;

          // std::cout<<"here it works 11"<<std::endl;

          const double aspect_timestep_in_years = this->get_timestep() / year_in_seconds;

          unsigned int fastscape_iterations = fastscape_steps_per_aspect_step;
          double fastscape_timestep_in_years = aspect_timestep_in_years / fastscape_iterations;
          while (fastscape_timestep_in_years > maximum_fastscape_timestep)
            {
              fastscape_iterations *= 2;
              fastscape_timestep_in_years *= 0.5;
            }

          std::cout << "[DEBUG] fastscape_steps_per_aspect_step = " << fastscape_steps_per_aspect_step << std::endl;
          std::cout << "[DEBUG] aspect_timestep_in_years = " << aspect_timestep_in_years << std::endl;
          std::cout << "[DEBUG] fastscape_iterations = " << fastscape_iterations << std::endl;
            

// TODO: This is a pretty arbitrary cell, with a pretty arbitrary cell diameter. Is that what you
// want? If so, perhaps say so in a comment and explain what it is used for.
//
// Separately, let's call this variable surface_cell_area in concordance with our convention
          const double surface_cell_area = surface_mesh.begin_active()->measure();

          // 1. Create the FastScape grid adapter
// TODO: I haven't looked at the GridAdapter yet, but we should see whether we can't get rid of the const_cast.
          auto grid = GridAdapterType(const_cast<SurfaceMeshType &>(surface_mesh), surface_cell_area);
          auto shape = grid.shape();  // Should be {n_grid_nodes}

          // 2. Define the node status array with correct shape and type
          xt::xtensor<fastscapelib::node_status, 1> node_status_array(shape);
          std::fill(node_status_array.begin(), node_status_array.end(), fastscapelib::node_status::core);
          grid.set_nodes_status(node_status_array);

          // // 3. Create the flow graph
          auto flow_graph = FlowGraphType(
              grid,
              {
                  std::make_shared<fastscapelib::single_flow_router>(),
                  std::make_shared<fastscapelib::mst_sink_resolver>()
              }
);

          // 4. Set base level nodes
          std::vector<std::size_t> base_level_nodes = {0};
          flow_graph.set_base_levels(base_level_nodes);

          // 5. Build the erosion model
          fastscapelib::spl_eroder<FlowGraphType> spl_eroder
          =
            // fastscapelib::make_spl_eroder(flow_graph, kff, n, m, 1e-5)
            fastscapelib::make_spl_eroder(flow_graph, kff, area_exp, slope_exp, 1e-5);
          std::cout << "KFF : " << kff<< std::endl;

          //Only for raster grid 
          //To propose if we add back the raster grid option later 
          // auto diffusion_eroder = fastscapelib::make_diffusion_adi_eroder(grid, kdd);

         // uplift rate in meters per year for Fastscape
          std::vector<double> uplift_rate_in_m_year(vz.size());
          for (size_t i = 0; i < vz.size(); ++i) {
              uplift_rate_in_m_year[i] = vz[i] / year_in_seconds;
          }

          // 6. Create data arrays using grid shape
          if (!elevation_initialized)
          {
            elevation = xt::adapt(h, shape);
            elevation_initialized = true;
            std::cout << "[DEBUG INIT] Initialized elevation from ASPECT solution.\n";
          }
          else
          {
            std::cout << "[DEBUG INIT] Reusing elevation from previous timestep.\n";
          }
          // auto elevation      = xt::adapt(h, shape);
          auto elevation_old  = xt::adapt(h_old, shape);
          auto uplift_rate    = xt::adapt(uplift_rate_in_m_year, shape);

          //To fix the borders for box models
          auto row_bounds = xt::view(uplift_rate, xt::keep(0, -1), xt::all());
          row_bounds = 0.0;
          auto col_bounds = xt::view(uplift_rate, xt::all(), xt::keep(0, -1));
          col_bounds = 0.0;

          xt::xarray<double> drainage_area  = xt::zeros<double>(shape);
          xt::xarray<double> sediment_flux  = xt::zeros<double>(shape);
          xt::xarray<double> spl_erosion = xt::zeros<double>(shape);
         xt::xarray<double> uplifted_elevation= xt::zeros<double>(shape);


          std::cout << "\n[DEBUG INIT] Max elevation = " << xt::amax(elevation)()
                    << ", Min = " << xt::amin(elevation)() << std::endl;
          std::cout << "[DEBUG INIT] Max uplift rate = " << xt::amax(uplift_rate)()
                    << ", Min = " << xt::amin(uplift_rate)() << std::endl;

    
// TODO: Is the comment wrong? We're not at the initial state any more, right?

// TODO: In the long run (perhaps not right away), we should adopt a system like many of
// the other plugins where we don't output information in *every* time step, but only
// when requested.


// --- FASTSCAPE LOOP ---
for (unsigned int i = 0; i < fastscape_iterations; ++i)
{
  std::cout << "\nFastScape iteration " << (i+1) << "/" << fastscape_iterations << std::endl;

  // Uplift
  uplifted_elevation = elevation + fastscape_timestep_in_years * uplift_rate;

  std::cout << "[DEBUG] uplifted_elevation max = " << xt::amax(uplifted_elevation)()
            << ", min = " << xt::amin(uplifted_elevation)() << std::endl;

  // Flow routing
    xt::noalias(drainage_area) = 0.0;

  flow_graph.update_routes(uplifted_elevation);
  flow_graph.accumulate(drainage_area, 1.0);

  std::cout << "[DEBUG] drainage_area max = " << xt::amax(drainage_area)()
            << ", min = " << xt::amin(drainage_area)() << std::endl;

  // Stream Power Law erosion
  auto spl_erosion = spl_eroder.erode(uplifted_elevation, drainage_area, fastscape_timestep_in_years);

  std::cout << "[DEBUG] spl_erosion max = " << xt::amax(spl_erosion)()
            << ", min = " << xt::amin(spl_erosion)() << std::endl;

  // Apply erosion
  sediment_flux = flow_graph.accumulate(spl_erosion);
  elevation = uplifted_elevation - spl_erosion;

  std::cout << "[DEBUG] updated elevation max = " << xt::amax(elevation)()
            << ", min = " << xt::amin(elevation)() << std::endl;
}


const unsigned int timestep_index = this->get_timestep_number();
const std::string output_directory = this->get_output_directory() + "/fastscapeCC";

if (!std::filesystem::exists(output_directory))
  std::filesystem::create_directories(output_directory);

// === Prepare data ===
dealii::Vector<double> elevation_output(n_grid_nodes);
dealii::Vector<double> uplift_rate_output(n_grid_nodes);
dealii::Vector<double> erosion_output(n_grid_nodes);
dealii::Vector<double> drainage_area_output(n_grid_nodes);

for (unsigned int j = 0; j < n_grid_nodes; ++j)
{
  elevation_output[j]     = elevation[j];
  uplift_rate_output[j]   = uplift_rate[j];
  erosion_output[j]       = spl_erosion[j];
  drainage_area_output[j] = drainage_area[j];
}

// === Write per-process VTU ===
dealii::DataOut<dim-1, dim> data_out;
data_out.attach_dof_handler(surface_mesh_dof_handler);
data_out.add_data_vector(elevation_output, "Elevation");
data_out.add_data_vector(uplift_rate_output, "UpliftRate");
data_out.add_data_vector(erosion_output, "Erosion");
data_out.add_data_vector(drainage_area_output, "DrainageArea");
data_out.build_patches();

const std::string basename = "fastscape_timestep_" + Utilities::int_to_string(timestep_index, 4);
const std::string filename_base = output_directory + "/" + basename;

const auto mpi_rank = Utilities::MPI::this_mpi_process(this->get_mpi_communicator());
const auto mpi_size = Utilities::MPI::n_mpi_processes(this->get_mpi_communicator());

const std::string vtu_filename = filename_base + "." + Utilities::int_to_string(mpi_rank, 4) + ".vtu";
{
  std::ofstream output(vtu_filename);
  data_out.write_vtu(output);
}

// === Rank 0: Write .pvtu and update .pvd ===
if (mpi_rank == 0)
{
  // .pvtu
  std::vector<std::string> vtu_filenames;
  for (unsigned int r = 0; r < mpi_size; ++r)
    vtu_filenames.push_back(basename + "." + Utilities::int_to_string(r, 4) + ".vtu");

  const std::string pvtu_filename = filename_base + ".pvtu";
  {
    std::ofstream pvtu(pvtu_filename);
    data_out.write_pvtu_record(pvtu, vtu_filenames);
  }

  // .pvd
  const std::string pvd_filename = output_directory + "/fastscape_timeseries.pvd";
  std::vector<std::string> existing_entries;
  if (std::filesystem::exists(pvd_filename))
  {
    std::ifstream in(pvd_filename);
    std::string line;
    while (std::getline(in, line))
      if (line.find("<DataSet") != std::string::npos)
        existing_entries.push_back(line);
  }

  std::ostringstream new_entry;
  new_entry << "    <DataSet timestep=\"" << timestep_index << "\" group=\"\" part=\"0\""
            << " file=\"" << basename << ".pvtu\"/>";
  existing_entries.push_back(new_entry.str());

  std::ofstream pvd_file(pvd_filename);
  pvd_file << "<?xml version=\"1.0\"?>\n"
           << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
           << "  <Collection>\n";
  for (const auto &entry : existing_entries)
    pvd_file << entry << "\n";
  pvd_file << "  </Collection>\n</VTKFile>\n";
}





          // // Create output directory
          // const std::string output_directory = (this->get_output_directory() + "/fastscapeCC");
          // if (!std::filesystem::exists(output_directory))
          //   std::filesystem::create_directories(output_directory);

          // // Vector to track output files
          // std::vector<std::string> vtu_filenames;


          // {
          //   dealii::Vector<double> elevation_output(n_grid_nodes);
          //   dealii::Vector<double> uplift_rate_output(n_grid_nodes);
          //   dealii::Vector<double> erosion_output(n_grid_nodes);
          //   dealii::Vector<double> drainage_area_output(n_grid_nodes);

          //   for (unsigned int j = 0; j < n_grid_nodes; ++j)
          //   {
          //     elevation_output[j]      = elevation[j];
          //     uplift_rate_output[j]    = uplift_rate[j];
          //     erosion_output[j]        = 0.0;
          //     drainage_area_output[j]  = 0.0;
          //   }

          //   dealii::DataOut<dim-1, dim> data_out;
          //   data_out.attach_dof_handler(surface_mesh_dof_handler);
          //   data_out.add_data_vector(elevation_output, "Elevation");
          //   data_out.add_data_vector(uplift_rate_output, "UpliftRate");
          //   data_out.add_data_vector(erosion_output, "Erosion");
          //   data_out.add_data_vector(drainage_area_output, "DrainageArea");
          //   data_out.build_patches();

          //   std::string filename = output_directory + "/fastscape_surface_iteration_0000.vtu";
          //   std::ofstream output(filename);
          //   data_out.write_vtu(output);
          //   this->get_pcout() << "➤ Wrote initial VTK file: " << filename << std::endl;

          //   vtu_filenames.push_back(filename);
          // }

          // // === FASTSCAPE LOOP ===
          // for (unsigned int i = 1; i <= fastscape_iterations; ++i)
          // {
          //   std::cout << "\nFastScape iteration " << i << "/" << fastscape_iterations << std::endl;

          //   xt::xarray<double> uplifted_elevation = elevation + fastscape_timestep_in_years * uplift_rate;

          //   std::cout << "[DEBUG] uplifted_elevation max = " << xt::amax(uplifted_elevation)()
          //             << ", min = " << xt::amin(uplifted_elevation)() << std::endl;

          //   flow_graph.update_routes(uplifted_elevation);
          //   flow_graph.accumulate(drainage_area, 1.0);

          //   std::cout << "[DEBUG] drainage_area max = " << xt::amax(drainage_area)()
          //             << ", min = " << xt::amin(drainage_area)() << std::endl;

          //   auto spl_erosion = spl_eroder.erode(uplifted_elevation, drainage_area, fastscape_timestep_in_years);

          //   std::cout << "[DEBUG] spl_erosion max = " << xt::amax(spl_erosion)()
          //             << ", min = " << xt::amin(spl_erosion)() << std::endl;

          //   //Working for raster grid only
          //   // auto diff_erosion = diffusion_eroder.erode(uplifted_elevation - spl_erosion, fastscape_timestep_in_years);

          //   sediment_flux = flow_graph.accumulate(spl_erosion);
          //   elevation = uplifted_elevation - spl_erosion ; //- diff_erosion for raster grid

          //   std::cout << "[DEBUG] updated elevation max = " << xt::amax(elevation)()
          //             << ", min = " << xt::amin(elevation)() << std::endl;

          //   // Prepare output vectors
          //   dealii::Vector<double> elevation_output(n_grid_nodes);
          //   dealii::Vector<double> uplift_rate_output(n_grid_nodes);
          //   dealii::Vector<double> erosion_output(n_grid_nodes);
          //   dealii::Vector<double> drainage_area_output(n_grid_nodes);

          //   for (unsigned int j = 0; j < n_grid_nodes; ++j)
          //   {
          //     elevation_output[j]      = elevation[j];
          //     uplift_rate_output[j]    = uplift_rate[j];
          //     erosion_output[j]        = spl_erosion[j];
          //     drainage_area_output[j]  = drainage_area[j];
          //   }

          //   dealii::DataOut<dim-1, dim> data_out;
          //   data_out.attach_dof_handler(surface_mesh_dof_handler);
          //   data_out.add_data_vector(elevation_output, "Elevation");
          //   data_out.add_data_vector(uplift_rate_output, "UpliftRate");
          //   data_out.add_data_vector(erosion_output, "Erosion");
          //   data_out.add_data_vector(drainage_area_output, "DrainageArea");
          //   data_out.build_patches();

          //   std::string filename = output_directory + "/fastscape_surface_iteration_" +
          //                         Utilities::int_to_string(i, 4) + ".vtu";
          //   std::ofstream output(filename);
          //   data_out.write_vtu(output);

          //   this->get_pcout() << "➤ Wrote VTK file: " << filename << std::endl;
          //   vtu_filenames.push_back(filename);
          // }

          // // === WRITE PVD FILE ===
          // {
          //   const std::string pvd_filename = output_directory + "/fastscape_iterations.pvd";
          //   std::ofstream pvd_file(pvd_filename);

          //   pvd_file << "<?xml version=\"1.0\"?>\n"
          //           << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
          //           << "  <Collection>\n";

          //   for (unsigned int i = 0; i < vtu_filenames.size(); ++i)
          //   {
          //     pvd_file << "    <DataSet timestep=\"" << i << "\" group=\"\" part=\"0\"\n"
          //             << "             file=\""
          //             << "fastscape_surface_iteration_" << Utilities::int_to_string(i, 4) << ".vtu\"/>\n";
          //   }

          //   pvd_file << "  </Collection>\n"
          //           << "</VTKFile>\n";

          //   this->get_pcout() << "➤ Wrote PVD file: " << pvd_filename << std::endl;
          // }
 
          // Compute erosion velocities
// TODO: We compute a vertical velocity in meters per year here. But internally, ASPECT
// computes in meters/second. When you apply these constraints, you'll have to multiple
// as appropriate. It might be nice to do that here already.
          const double aspect_timestep_in_seconds = aspect_timestep_in_years * year_in_seconds;

          for (unsigned int i = 0; i < n_grid_nodes; ++i)
            {
              V[i] = (elevation[i] - elevation_old[i]) / aspect_timestep_in_seconds;
              // std::cout << "grid node: " << i << ", V: " << V[i] << std::endl;
            }


          std::cout<<"here it works 13"<<std::endl;

          // Broadcast V to all other processes
          MPI_Bcast(&V[0], n_grid_nodes, MPI_DOUBLE, 0, this->get_mpi_communicator());
        }
      else
        {
// TODO: This is the else branch I mention above
          for (unsigned int i = 0; i < temporary_variables.size(); ++i)
            MPI_Ssend(&temporary_variables[i][0], temporary_variables[1].size(), MPI_DOUBLE, 0, 42, this->get_mpi_communicator());

// TODO: This one corresponds to the Bcast at the end of the if-block above and sends all other processes
// what process 0 has computed. I would move that *after* the if-block in the same
// way as the two lines above were moved *before* the if-block, as that more clearly
// communicates the order of communication.
          MPI_Bcast(&V[0], n_grid_nodes, MPI_DOUBLE, 0, this->get_mpi_communicator());
        }

      // Maybe this is enough
      // Step 1: Get V ordered
      auto erosion_function = [&](const Point<dim> &p) -> double
      {
        const unsigned int index = this->vertex_index(p);
        return V[index]; //meters/year for ASPECT
      };
      std::cout << "here it works 14" << std::endl;

      //Eventually we could interpolate for better accuracy
      // auto erosion_function = [&](const Point<dim> &p) -> double
      // {
      //   try
      //   {
      //     const auto cell_and_point = dealii::GridTools::find_active_cell_around_point<dim - 1>(*surface_cache, p);

      //     const auto &cell = cell_and_point.first;

      //     if (cell->is_valid())
      //     {
      //       // Interpolate using FEValues
      //       const MappingQ1<dim - 1> mapping;
      //       FEValues<dim - 1> fe_values(mapping,
      //                                   surface_fe,
      //                                   Quadrature<dim - 1>(cell_and_point.second),
      //                                   update_values);

      //       fe_values.reinit(cell);

      //       std::vector<double> values(fe_values.dofs_per_cell);
      //       cell->get_dof_values(V, values); 

      //       return values[0];
      //     }
      //     else
      //     {
      //       return 0.0;
      //     }
      //   }
      //   catch (...)
      //   {
      //     // Fallback for out-of-bound or failure
      //     return 0.0;
      //   }
      // };


      VectorFunctionFromScalarFunctionObject<dim> radial_velocity_field(
        erosion_function,
        dim - 1, // project onto radial component
        dim      // full space dimension
      );
      std::cout << "Rank " << Utilities::MPI::this_mpi_process(this->get_mpi_communicator())
                << ": here it works 15" << std::endl;



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
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Box");
        {
          prm.declare_entry ("X repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in X direction.");
          prm.declare_entry ("Y repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in Y direction.");
          prm.declare_entry ("Z repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in Z direction.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

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
          prm.declare_entry("Surface refinement level", "1",
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

          prm.declare_entry("Output internal fastscape steps", "false",
                  Patterns::Bool(),
                  "Whether to output visualization files for internal FastScape steps.");


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
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Box");
        {
          repetitions[0] = prm.get_integer ("X repetitions");
          if (dim >= 2)
            {
              repetitions[1] = prm.get_integer ("Y repetitions");
            }
          if (dim >= 3)
            {
              repetitions[dim-1] = prm.get_integer ("Z repetitions");
            }
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection("Fastscapecc");
        {
          fastscape_steps_per_aspect_step = prm.get_integer("Number of steps");
          maximum_fastscape_timestep = prm.get_double("Maximum timestep");
          additional_refinement_levels = prm.get_integer("Additional fastscape refinement levels");
          surface_refinement_level = prm.get_integer("Surface refinement level");
          surface_refinement_difference = prm.get_integer("Surface refinement difference");
          noise_h = prm.get_double("Initial noise magnitude");
          output_internal_fastscape_steps = prm.get_bool("Output internal fastscape steps");


          if (!this->convert_output_to_years())
            {
              maximum_fastscape_timestep /= year_in_seconds;
            }

          prm.enter_subsection("Erosional parameters");
          {
            area_exp = prm.get_double("Drainage area exponent");
            slope_exp = prm.get_double("Slope exponent");
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

#endif // #ifdef ASPECT_WITH_FASTSCAPELIB
