/*
  Copyright (C) 2018 - 2023 by the authors of the ASPECT code.

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



#include <deal.II/numerics/vector_tools.h>
#include <aspect/mesh_deformation/openlem.h>
#include <aspect/geometry_model/box.h>
namespace aspect
{
  namespace MeshDeformation
  {
    template <int dim>
    OpenLEM<dim>::OpenLEM()
    {}


    template <int dim>
    void
    OpenLEM<dim>::initialize ()
    {
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          CitationInfo::add("openlem");

          AssertThrow(Plugins::plugin_type_matches<const GeometryModel::Box<dim>>(this->get_geometry_model()),
                      ExcMessage("OpenLEM can only be run with a box geometry model."));

          const GeometryModel::Box<dim> *geometry
            = dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model());

          // Find the id associated with the top boundary and boundaries that call mesh deformation.
          const types::boundary_id top_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id ("top");
          const std::set<types::boundary_id> mesh_deformation_boundary_ids
            = this->get_mesh_deformation_handler().get_active_mesh_deformation_boundary_indicators();

          // Get the deformation type names called for each boundary.
          std::map<types::boundary_id, std::vector<std::string>> mesh_deformation_boundary_indicators_map
            = this->get_mesh_deformation_handler().get_active_mesh_deformation_names();

          // Loop over each mesh deformation boundary, and make sure openlem is only called on the surface.
          for (const types::boundary_id id : mesh_deformation_boundary_ids)
            {
              const std::vector<std::string> &names = mesh_deformation_boundary_indicators_map[id];
              for (const auto &name : names)
                {
                  if (name == "openlem")
                    AssertThrow(id == top_boundary,
                                ExcMessage("OpenLEM can only be called on the surface boundary."));
                }
            }

          // Several compositional fields are commonly used in conjunction with the openlem plugin, i.e.
          // "sediment_age" to track the age of the sediment deposited and "deposition_depth" to track the depth
          // with respect to the unperturbed surface of the model domain. Their values are controlled by setting
          // boundary conditions on the top boundary that is deformed by openlem. While it is useful to track these
          // fields, they are not needed for any function in the openlem plugin. If they exist however, we need
          // to make sure that these fields do not have the type "chemical composition" and are therefore not taken
          // into account when computing material properties.
          const std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();
          if (this->introspection().compositional_name_exists("sediment_age"))
            {
              const std::vector<std::string>::const_iterator
              it = std::find(chemical_field_names.begin(), chemical_field_names.end(), "sediment_age");
              AssertThrow (it == chemical_field_names.end(),
                           ExcMessage("There is a field sediment_age that is of type chemical composition. "
                                      "Please change it to type generic so that it does not affect material properties."));
            }
          if (this->introspection().compositional_name_exists("deposition_depth"))
            {
              const std::vector<std::string>::const_iterator
              it = std::find(chemical_field_names.begin(), chemical_field_names.end(), "deposition_depth");
              AssertThrow (it == chemical_field_names.end(),
                           ExcMessage("There is a field deposition_depth that is of type chemical composition. "
                                      "Please change it to type generic so that it does not affect material properties."));
            }

          // Initialize parameters for restarting openlem
          restart = this->get_parameters().resume_computation;

          // Since we don't open these until we're on one process, we need to check if the
          // restart files exist beforehand.
          if (restart)
            {
              if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
                {
                  AssertThrow(Utilities::fexists(this->get_output_directory() + "openlem_elevation_restart.txt"),
                              ExcMessage("Cannot open topography file to restart openLEM."));
                  AssertThrow(Utilities::fexists(this->get_output_directory() + "openlem_basement_restart.txt"),
                              ExcMessage("Cannot open basement file to restart openLEM."));
                  AssertThrow(Utilities::fexists(this->get_output_directory() + "openLEM_silt_fraction_restart.txt"),
                              ExcMessage("Cannot open silt fraction file to restart openLEM."));
                }
            }
          // The first entry represents the minimum coordinates of the model domain, the second the model extent.
          for (unsigned int d=0; d<dim; ++d)
            {
              grid_extent[d].first = geometry->get_origin()[d];
              grid_extent[d].second = geometry->get_extents()[d];
            }

          // Get the x and y repetitions used in the parameter file so
          // the openlem cell size can be properly set.
          const std::array<unsigned int, dim> repetitions = geometry->get_repetitions();

          // Set number of x points, which is generally 1+(openlem refinement level)^2.
          // The openlem refinement level is a combination of the maximum ASPECT refinement level
          // at the surface and any additional refinement we want in openlem. If
          // repetitions are specified we need to adjust the number of points to match what ASPECT has,
          // which can be determined by multiplying the points by the repetitions before adding 1.
          // Finally, if ghost nodes are used we add two additional points on each side.
          const unsigned int ghost_nodes = 2*use_ghost_nodes;
          const unsigned int openlem_refinement_level = maximum_surface_refinement_level + additional_refinement_levels;
          const unsigned int openlem_nodes = Utilities::pow(2,openlem_refinement_level);
          openlem_nx = openlem_nodes * repetitions[0] + ghost_nodes + 1;

          // Size of openlem cell.
          openlem_dx = (grid_extent[0].second)/(openlem_nodes * repetitions[0]);

          // openlem X extent, which is generally ASPECT's extent unless the ghost nodes are used,
          // in which case 2 cells are added on either side.
          openlem_x_extent = (grid_extent[0].second) + openlem_dx * ghost_nodes;

          // Sub intervals are 3 less than points, if including the ghost nodes. Otherwise 1 less.
          table_intervals[0] = openlem_nodes * repetitions[0];
          table_intervals[dim-1] = 1;

          if (dim == 2)
            {
              openlem_dy = openlem_dx;
              openlem_y_extent = std::round(openlem_y_extent_2d/openlem_dy)*openlem_dy + openlem_dy * ghost_nodes;
              openlem_ny = 1+openlem_y_extent/openlem_dy;
            }
          else
            {
              openlem_ny = openlem_nodes * repetitions[1] + ghost_nodes + 1;
              openlem_dy = (grid_extent[1].second)/(openlem_nodes * repetitions[1]);
              table_intervals[1] = openlem_nodes * repetitions[1];
              openlem_y_extent = (grid_extent[1].second) + openlem_dy * ghost_nodes;
            }


          // Create a folder for the openlem visualization files.
          Utilities::create_directory (this->get_output_directory() + "openLEM/",
                                       this->get_mpi_communicator(),
                                       false);

          last_output_time = 0;
  	  // preparte the openLEM variables
          grid_old = openlem::Grid<>(openlem_nx,openlem_ny);
          grid_new = openlem::Grid<>(openlem_nx,openlem_ny);

	  // todo: fill the grid new

          // Define all nodes with non-positive elevations (here, the ocean around the
          // island as boundary nodes
          for ( int i = 0; i < grid_new.m; ++i )
            for ( int j = 0; j < grid_new.n; ++j )
              grid_new[i][j].b = grid_new[i][j].h <= 0;

	  // Compute initial flow pattern and water level
	  grid_new.fillLakes();

        }
    }

    template <int dim>
    void
    OpenLEM<dim>::update ()
    {

      // Because there is no increase in time during timestep 0, we return and only
      // initialize and run openlem from timestep 1 and on.
      if (this->get_timestep_number() == 0)
        return;

      TimerOutput::Scope timer_section(this->get_computing_timer(), "openlem plugin");

      const unsigned int current_timestep = this->get_timestep_number ();
      const double aspect_timestep_in_years = this->get_timestep() / year_in_seconds;

      // Find a openlem timestep that is below our maximum timestep.
      unsigned int openlem_iterations = openlem_steps_per_aspect_step;
      double openlem_timestep_in_years = aspect_timestep_in_years/openlem_iterations;
      while (openlem_timestep_in_years>maximum_openlem_timestep)
        {
          openlem_iterations *= 2;
          openlem_timestep_in_years *= 0.5;
        }

      // Vector to hold the velocities that represent the change to the surface.
      const unsigned int openlem_array_size = openlem_nx*openlem_ny;
      std::vector<double> mesh_velocity_z(openlem_array_size);

      // openlem requires multiple specially defined and ordered variables sent to its functions. To make
      // the transfer of these down to one process easier, we first fill out a vector of local_aspect_values,
      // then when we get down to one process we use these local_aspect_values to fill the double arrays
      // in the order needed for openlem.
      std::vector<std::vector<double>> local_aspect_values = get_aspect_values();

      // Run openlem on single process.
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          // Initialize the variables that will be sent to openlem.
          // Elevation is initialized at a very high number so that we can later check that all points
          // received data from ASPECT, and if not throw an assert.

          fill_openlem_arrays(grid_new,
                                 local_aspect_values);

          if (current_timestep == 1 || restart)
            {
              this->get_pcout() << "   Initializing openlem... " << (1+maximum_surface_refinement_level+additional_refinement_levels) <<
                                " levels, cell size: " << openlem_dx << " m." << std::endl;


              //initialize_openlem(elevation,
              //                     basement,
              //                     bedrock_transport_coefficient_array,
              //                     bedrock_river_incision_rate_array,
              //                     silt_fraction);
            }
          else
            {
              // If it isn't the first timestep we ignore initialization and instead copy all height values from openlem.
              // Generally, we overwrite the topography data from ASPECT as openlem may be at a higher resolution. However,
              // if we are not using openlem to advect then we do not want to do this and instead use the ASPECT values.
              // if (openlem_advection_uplift)
              //   openlem_copy_h_(elevation.data());
            }

          // Find the appropriate sediment rain based off the time interval.
          const double time_in_years = this->get_time() / year_in_seconds;
          //auto it = std::lower_bound(sediment_rain_times.begin(), sediment_rain_times.end(), time_in_years);
          //const unsigned int inds = std::distance(sediment_rain_times.begin(), it);
          //const double sediment_rain = sediment_rain_rates[inds];

          // Keep initial h values so we can calculate velocity later.
          // In the first timestep, h will be given from other processes.
          // In later timesteps, we copy h directly from openlem.
          //std::mt19937 random_number_generator(openlem_seed);
          //std::uniform_real_distribution<double> random_distribution(-noise_elevation,noise_elevation);
          for (unsigned int i=0; i<openlem_array_size; ++i)
            {
              //elevation_old[i] = elevation[i];

              // Initialize random noise after elevation_old is set, so ASPECT sees this initial topography change.
              // Changing boundary height directly on a fixed openlem boundary causes reproducibility issues,
              // as such we do not add noise to the boundaries regardless of whether they are ghost nodes
              // or not. However, the boundaries can be changed using the uplift velocity and not cause
              // these issues.
              // TODO: Should this be done through velocities instead of a flat height change?
              // if (!is_ghost_node(i,true))
              //   {
              //     if (current_timestep == 1)
              //       {
              //         // + or - topography based on the initial noise magnitude.
              //         //const double elevation_seed = random_distribution(random_number_generator);
              //         //elevation[i] = elevation[i] + elevation_seed;
              //       }

              //     // Here we add the sediment rain (m/yr) as a flat increase in height.
              //     // This is done because adding it as an uplift rate would affect the basement.
              //     // if (sediment_rain > 0 && use_marine_component)
              //     //   {
              //     //     // Only apply sediment rain to areas below sea level.
              //     //     if (elevation[i] < sea_level)
              //     //       {
              //     //         // If the rain would put us above sea level, set height to sea level.
              //     //         if (elevation[i] + sediment_rain*aspect_timestep_in_years > sea_level)
              //     //           elevation[i] = sea_level;
              //     //         else
              //     //           elevation[i] = std::min(sea_level,elevation[i] + sediment_rain*aspect_timestep_in_years);
              //     //       }
              //     //   }
              //   }
            }

          // The ghost nodes are added as a single layer of points surrounding the entire model.
          // For example, if ASPECT's surface mesh is a 2D surface that is 3x3 (nx x ny) points,
          // openlem will be set as a 2D 5x5 point surface. On return to ASPECT, the outer ghost nodes
          // will be ignored, and ASPECT will see only the inner 3x3 surface of openlem.
          //if (use_ghost_nodes)
          //  set_ghost_nodes(elevation,
          //                  velocity_x,
          //                  velocity_y,
          //                  velocity_z,
          //                  openlem_timestep_in_years,
          //                  false);

          // If specified, apply the orographic controls to the openlem model.
          // if (use_orographic_controls)
          //   apply_orographic_controls(elevation,
          //                             bedrock_transport_coefficient_array,
          //                             bedrock_river_incision_rate_array);

          // Set velocity components.
          //if (openlem_advection_uplift)
          //  {
          //    //openlem_set_u_(velocity_z.data());
          //    //openlem_set_v_(velocity_x.data(),
          //                     velocity_y.data());
          //  }

          // Set h to new values, and erosional parameters if there have been changes.
          //openlem_set_h_(elevation.data());

          // openlem_set_erosional_parameters_(bedrock_river_incision_rate_array.data(),
          //                                     &sediment_river_incision_rate,
          //                                     &drainage_area_exponent_m,
          //                                     &slope_exponent_n,
          //                                     bedrock_transport_coefficient_array.data(),
          //                                     &sediment_transport_coefficient,
          //                                     &bedrock_deposition_g,
          //                                     &sediment_deposition_g,
          //                                     &slope_exponent_p);

          // Find timestep size, run openlem, and make visualizations.
          //execute_openlem(elevation,
          //                  bedrock_transport_coefficient_array,
          //                  velocity_x,
          //                  velocity_y,
          //                  velocity_z,
          //                  openlem_timestep_in_years,
          //                  openlem_iterations);

          // Write a file to store h, b & step for restarting.
          // TODO: It would be good to roll this into the general ASPECT checkpointing,
          // and when we do this needs to be changed.
          if (((this->get_parameters().checkpoint_time_secs == 0) &&
               (this->get_parameters().checkpoint_steps > 0) &&
               ((current_timestep + 1) % this->get_parameters().checkpoint_steps == 0)) ||
              (this->get_time() == this->get_end_time() && this->get_timestepping_manager().need_checkpoint_on_terminate()))
            {
              //save_restart_files(elevation,
              //                   basement,
              //                   silt_fraction);
            }

          // Find out our velocities from the change in height.
          // Where mesh_velocity_z is a vector of array size that exists on all processes.
          for (unsigned int ix=0; ix<openlem_nx; ++ix)
          for (unsigned int iy=0; iy<openlem_ny; ++iy)
            {
	      openlem::Node* node = grid.getNode(index_x,index_y);
              mesh_velocity_z[i] = (elevation[i] - elevation_old[i])/aspect_timestep_in_years;
            }

          Utilities::MPI::broadcast(this->get_mpi_communicator(), mesh_velocity_z, 0);
        }
      else
        {
          for (unsigned int i=0; i<local_aspect_values.size(); ++i)
            MPI_Ssend(&local_aspect_values[i][0], local_aspect_values[1].size(), MPI_DOUBLE, 0, 42, this->get_mpi_communicator());

          // Check whether the openlem mesh was filled with data.
          const bool openlem_mesh_filled = Utilities::MPI::broadcast (this->get_mpi_communicator(), true, 0);
          if (openlem_mesh_filled != true)
            throw aspect::QuietException();

          // This is called solely so we can set the timer and will return immediately.
          execute_openlem(mesh_velocity_z,
                            mesh_velocity_z,
                            mesh_velocity_z,
                            mesh_velocity_z,
                            mesh_velocity_z,
                            aspect_timestep_in_years,
                            openlem_steps_per_aspect_step);

          mesh_velocity_z = Utilities::MPI::broadcast(this->get_mpi_communicator(), mesh_velocity_z, 0);
        }
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years

      // update
      // TODO: get velocity at surface from all processes.
      const double dt = 1e-3; //year?
      for(int i = 0; i < 100; ++i )
      {
        int nc = grid_new.computeFlowDirection();
        printf("Changes in flow direction: %i\n",nc);  
        double ch = grid_new.erode(dt);
        printf("Maximum elevation change: %e\n",ch);
      }
       // Todo: height_diff = grid_new->height - grid_old->height;
       // Todo: interpolate hight_diff to velocties on nodes
       // Todo: move info to all other processes
       // Todo: copy grid_new to grid_old
    }



    /**
     * A function that creates constraints for the velocity of certain mesh
     * vertices (e.g. the surface vertices) for a specific boundary.
     * The calling class will respect
     * these constraints when computing the new vertex positions.
     */
    template <int dim>
    void
    OpenLEM<dim>::compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                           AffineConstraints<double> &mesh_velocity_constraints,
                                                           const std::set<types::boundary_id> &boundary_ids) const
    {
      // Loop over all boundary indicators to set the velocity constraints
      /*for (const auto boundary_id : boundary_ids)
        VectorTools::interpolate_boundary_values (this->get_mapping(),
                                                  mesh_deformation_dof_handler,
                                                  boundary_id,
                                                  function,
                                                  mesh_velocity_constraints);
      */  

      // Get the sizes needed for a data table of the mesh velocities.
      TableIndices<dim> size_idx;
      for (unsigned int d=0; d<dim; ++d)
        {
          size_idx[d] = table_intervals[d]+1;
        }

      // Initialize a table to hold all velocity values that will be interpolated back to ASPECT.
      const Table<dim,double> velocity_table = fill_data_table(mesh_velocity_z, size_idx, openlem_nx, openlem_ny);

      // As our grid_extent variable end points do not account for the change related to an origin
      // not at 0, we adjust this here into an interpolation extent.
      std::array<std::pair<double,double>,dim> interpolation_extent;
      for (unsigned int d=0; d<dim; ++d)
        {
          interpolation_extent[d].first = grid_extent[d].first;
          interpolation_extent[d].second = (grid_extent[d].second + grid_extent[d].first);
        }

      //Functions::InterpolatedUniformGridData<dim> *velocities;
      Functions::InterpolatedUniformGridData<dim> velocities (interpolation_extent,
                                                              table_intervals,
                                                              velocity_table);

      VectorFunctionFromScalarFunctionObject<dim> vector_function_object(
        [&](const Point<dim> &p) -> double
      {
        return velocities.value(p);
      },
      dim-1,
      dim);

      VectorTools::interpolate_boundary_values (mesh_deformation_dof_handler,
                                                *boundary_ids.begin(),
                                                vector_function_object,
                                                mesh_velocity_constraints);
    
    }

    template <int dim>
    std::vector<std::vector<double>>
    OpenLEM<dim>::get_aspect_values() const 
{

      const types::boundary_id relevant_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id ("top");
      std::vector<std::vector<double>> local_aspect_values(dim+2, std::vector<double>());

      // Get a quadrature rule that exists only on the corners, and increase the refinement if specified.
      const QIterated<dim-1> face_corners (QTrapezoid<1>(),
                                           Utilities::pow(2,additional_refinement_levels+surface_refinement_difference));

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        face_corners,
                                        update_values |
                                        update_quadrature_points);

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned() && cell->at_boundary())
          for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
            if (cell->face(face_no)->at_boundary())
              {
                if ( cell->face(face_no)->boundary_id() != relevant_boundary)
                  continue;

                std::vector<Tensor<1,dim>> vel(face_corners.size());
                fe_face_values.reinit(cell, face_no);
                fe_face_values[this->introspection().extractors.velocities].get_function_values(this->get_solution(), vel);

                for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
                  {
                    const Point<dim> vertex = fe_face_values.quadrature_point(corner);

                    // Find what x point we're at. Add 1 or 2 depending on if ghost nodes are used.
                    // Subtract the origin point so that it corresponds to an origin of 0,0 in openlem.
                    const double indx = (vertex(0) - grid_extent[0].first)/openlem_dx;

                    // The quadrature rule is created so that there are enough interpolation points in the
                    // lowest resolved ASPECT surface cell to fill out the openlem mesh. However, as the
                    // same rule is used for all cell sizes, higher resolution areas will have interpolation
                    // points that do not correspond to a openlem node. In which case, indx will not be a
                    // whole number and we can ignore the point.
                    if (std::abs(indx - std::round(indx)) >= node_tolerance)
                      continue;


                    // If we're in 2D, we want to take the values and apply them to every row of X points.
                    if (dim == 2)
                      {
                        for (unsigned int ys=0; ys<openlem_ny; ++ys)
                          {
                            // openlem indexes from 1 to n, starting at X and Y = 0, and increases
                            // across the X row. At the end of the row, it jumps back to X = 0
                            // and up to the next X row in increasing Y direction. We track
                            // this to correctly place the variables later on.
                            // Nx*ys effectively tells us what row we are in
                            // and then indx tells us what position in that row.
                            //const double index = std::round(indx)+openlem_nx*ys;

                            local_aspect_values[0].push_back(vertex(dim-1) - grid_extent[dim-1].second);
                            local_aspect_values[1].push_back(indx);
                            local_aspect_values[2].push_back(ys);

                            for (unsigned int d=0; d<dim; ++d)
                              {
                                // Always convert to m/yr for openlem
                                local_aspect_values[3+d].push_back(vel[corner][d]*year_in_seconds);
                              }
                          }
                      }
                    // 3D case
                    else
                      {
                        // Because indy only gives us the row we're in, we don't need to add 2 for the ghost node.
                        const double indy = 1+use_ghost_nodes+(vertex(1) - grid_extent[1].first)/openlem_dy;

                        if (std::abs(indy - std::round(indy)) >= node_tolerance)
                          continue;

                        //const double index = std::round((indy-1))*openlem_nx+std::round(indx);

                        local_aspect_values[0].push_back(vertex(dim-1) - grid_extent[dim-1].second);   //z component
                        local_aspect_values[1].push_back(indx);
                        local_aspect_values[2].push_back(indy);

                        for (unsigned int d=0; d<dim; ++d)
                          {
                            local_aspect_values[3+d].push_back(vel[corner][d]*year_in_seconds);
                          }
                      }
                  }
              }

      return local_aspect_values;
    }
    template <int dim>
    void OpenLEM<dim>::execute_openlem(std::vector<double> &elevation,
                                           std::vector<double> &extra_vtk_field,
                                           std::vector<double> &velocity_x,
                                           std::vector<double> &velocity_y,
                                           std::vector<double> &velocity_z,
                                           const double &openlem_timestep_in_years,
                                           const unsigned int &openlem_iterations) const
{
      TimerOutput::Scope timer_section(this->get_computing_timer(), "Execute openlem");
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) != 0)
        return;

      // Because on the first timestep we will create an initial VTK file before running openlem
      // and a second after, we first set the visualization step to zero.
      unsigned int visualization_step = 0;
      const unsigned int current_timestep = this->get_timestep_number ();
      std::string dirname = (this->get_output_directory() + "openlem/");
      const char *dirname_char=dirname.c_str();
      const unsigned int dirname_length = dirname.length();

      // Set time step
      //openlem_set_dt_(&openlem_timestep_in_years);
      this->get_pcout() << "   Executing openlem... " << (openlem_iterations) << " timesteps of " << openlem_timestep_in_years << " years." << std::endl;
      {
        // If it is the first timestep, write an initial VTK file.
        if (current_timestep == 1)
          {
            this->get_pcout() << "      Writing initial VTK..." << std::endl;
            // openlem by default visualizes a field called HHHHH,
            // and the parameter this shows will be whatever is given as the first
            // position. At the moment it visualizes the bedrock diffusivity.
            //openlem_named_vtk_(extra_vtk_field.data(),
            //                     &vexp,
            //                     &visualization_step,
            //                     dirname_char,
            //                     &dirname_length);
          }

        for (unsigned int openlem_iteration = 0; openlem_iteration < openlem_iterations; ++openlem_iteration)
          {
            //openlem_execute_step_();

            //// If we are using the ghost nodes we want to reset them every openlem timestep.
            //if (use_ghost_nodes)
            //  {
            //    openlem_copy_h_(elevation.data());

            //    set_ghost_nodes(elevation,
            //                    velocity_x,
            //                    velocity_y,
            //                    velocity_z,
            //                    openlem_timestep_in_years,
            //                    false);

            //    // Set velocity components.
            //    if (openlem_advection_uplift)
            //      {
            //        openlem_set_u_(velocity_z.data());
            //        openlem_set_v_(velocity_x.data(),
            //                         velocity_y.data());
            //      }

            //    // Set h to new values, and erosional parameters if there have been changes.
            //    openlem_set_h_(elevation.data());
            //  }
          }

        // Copy h values.
        //openlem_copy_h_(elevation.data());


        // Determine whether to create a VTK file this timestep.
        bool write_vtk = false;

        // if (this->get_time() >= last_output_time + output_interval || this->get_time() == this->get_end_time())
        //   {
        //     write_vtk = true;

        //     if (output_interval > 0)
        //       {
        //         // We need to find the last time output was supposed to be written.
        //         // this is the last_output_time plus the largest positive multiple
        //         // of output_intervals that passed since then. We need to handle the
        //         // edge case where last_output_time+output_interval==current_time,
        //         // we did an output and std::floor sadly rounds to zero. This is done
        //         // by forcing std::floor to round 1.0-eps to 1.0.
        //         const double magic = 1.0+2.0*std::numeric_limits<double>::epsilon();
        //         last_output_time = last_output_time + std::floor((this->get_time()-last_output_time)/output_interval*magic) * output_interval/magic;
        //       }
        //   }

        // if (write_vtk)
        //   {
        //     this->get_pcout() << "      Writing openlem VTK..." << std::endl;
        //     visualization_step = current_timestep;
        //     openlem_named_vtk_(extra_vtk_field.data(),
        //                          &vexp,
        //                          &visualization_step,
        //                          dirname_char,
        //                          &dirname_length);
        //   }
      }
    }


    template <int dim>
    void OpenLEM<dim>::fill_openlem_arrays(openlem::Grid<openlem::Node> &grid,
					       // std::vector<double> &elevation,
                                               // std::vector<double> &bedrock_transport_coefficient_array,
                                               // std::vector<double> &bedrock_river_incision_rate_array,
                                               // std::vector<double> &velocity_x,
                                               // std::vector<double> &velocity_y,
                                               // std::vector<double> &velocity_z,
                                               std::vector<std::vector<double>> &local_aspect_values) 
{
      for (unsigned int i=0; i<local_aspect_values[1].size(); ++i)
        {

          int index_x = local_aspect_values[1][i];
          int index_y = local_aspect_values[2][i];
	  openlem::Node* node = grid.getNode(index_x,index_y);
	  node->h = local_aspect_values[0][i];
	  node->l = 0;
	  node->u = local_aspect_values[dim+2][i];
          //elevation[index] = local_aspect_values[0][i];
          //velocity_x[index] = local_aspect_values[2][i];
          //velocity_z[index] = local_aspect_values[dim+1][i];

          //if (dim == 2)
          //  velocity_y[index] = 0;
          //else
          //  velocity_y[index] = local_aspect_values[3][i];
        }

      for (unsigned int p=1; p<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++p)
        {
          // First, find out the size of the array a process wants to send.
          MPI_Status status;
          MPI_Probe(p, 42, this->get_mpi_communicator(), &status);
          int incoming_size = 0;
          MPI_Get_count(&status, MPI_DOUBLE, &incoming_size);

          // Resize the array so it fits whatever the process sends.
          for (unsigned int i=0; i<local_aspect_values.size(); ++i)
            {
              local_aspect_values[i].resize(incoming_size);
            }

          for (unsigned int i=0; i<local_aspect_values.size(); ++i)
            MPI_Recv(&local_aspect_values[i][0], incoming_size, MPI_DOUBLE, p, 42, this->get_mpi_communicator(), &status);

          // Now, place the numbers into the correct place based off the index.
          for (unsigned int i=0; i<local_aspect_values[1].size(); ++i)
            {
              int index_x = local_aspect_values[1][i];
              int index_y = local_aspect_values[2][i];
	      openlem::Node* node = grid.getNode(index_x,index_y);
	      node->h = local_aspect_values[0][i];
	      node->l = 0;
	      node->u = local_aspect_values[dim+2][i];
              //const unsigned int index = local_aspect_values[1][i];
              //elevation[index] = local_aspect_values[0][i];
              //velocity_x[index] = local_aspect_values[2][i];
              //velocity_z[index] = local_aspect_values[dim+1][i];

              //// In 2D there are no y velocities, so we set them to zero.
              //if (dim == 2 )
              //  velocity_y[index] = 0;
              //else
              //  velocity_y[index] = local_aspect_values[3][i];
            }
        }

      // Initialize the bedrock river incision rate and transport coefficient,
      // and check that there are no empty mesh points due to
      // an improperly set maximum_surface_refinement_level, additional_refinement_levels,
      // and surface_refinement_difference
      bool openlem_mesh_filled = true;
      const unsigned int openlem_array_size = openlem_nx*openlem_ny;
      for (unsigned int i=0; i<openlem_array_size; ++i)
        {
          //bedrock_river_incision_rate_array[i] = bedrock_river_incision_rate;
          //bedrock_transport_coefficient_array[i] = bedrock_transport_coefficient;

          //// If this is a boundary node that is a ghost node then ignore that it
          //// has not filled yet as the ghost nodes haven't been set.
          //if (elevation[i] == std::numeric_limits<double>::max() && !is_ghost_node(i,false))
          //  openlem_mesh_filled = false;
        }

      Utilities::MPI::broadcast(this->get_mpi_communicator(), openlem_mesh_filled, 0);
      AssertThrow (openlem_mesh_filled == true,
                   ExcMessage("The openlem mesh is missing data. A likely cause for this is that the "
                              "maximum surface refinement or surface refinement difference are improperly set."));
    }



    template <int dim>
    bool
    OpenLEM<dim>::
    needs_surface_stabilization () const
    {
      return false;
    }



    template <int dim>
    void OpenLEM<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection ("openLEM");
        {
        }
        prm.leave_subsection();
      }
      prm.leave_subsection ();
    }

    template <int dim>
    void OpenLEM<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection("openLEM");
        {
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
    ASPECT_REGISTER_MESH_DEFORMATION_MODEL(OpenLEM,
                                           "openLEM",
                                           "TODO: A plugin, which prescribes the surface mesh to "
                                           "deform according to an analytically prescribed "
                                           "function. Note that the function prescribes a "
                                           "deformation velocity, i.e. the return value of "
                                           "this plugin is later multiplied by the time step length "
                                           "to compute the displacement increment in this time step. "
                                           "Although the function's time variable is interpreted as "
                                           "years when Use years in output instead of seconds is set to true, "
                                           "the boundary deformation velocity should still be given "
                                           "in m/s. The format of the "
                                           "functions follows the syntax understood by the "
                                           "muparser library, see {ref}\\`sec:run-aspect:parameters-overview:muparser-format\\`.")
  }
}
