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
#include <aspect/global.h>

#ifdef ASPECT_WITH_FASTSCAPE

#include <aspect/mesh_deformation/fastscape.h>
#include <aspect/geometry_model/box.h>
#include <deal.II/numerics/vector_tools.h>
#include <aspect/postprocess/visualization.h>
#include <ctime>
#include <aspect/simulator.h>

namespace aspect
{
  namespace MeshDeformation
  {
    /**
     * Define FastScape functions as C functions. Must use the exact same function/variable name
     * and type as used in FastScape. All function names must be made lowercase, and an
     * underscore added at the end. Types must be defined as pointers, and sent to
     * FastScape as a reference. Additional functions are available within FastScape,
     * see https://fastscape.org/fastscapelib-fortran/ for a list of all functions and
     * their input parameters. These functions must be defined at the top here before
     * they are used.
     */
    extern"C"
    {
      /**
       * Function to initialize FastScape.
       */
      void fastscape_init_();

      /**
       * Set the x and y extent of the FastScape model.
       */
      void fastscape_set_xl_yl_(const double *xxl,
                                const double *yyl);

      /**
       * Set number of grid points in x (nx) and y (ny)
       */
      void fastscape_set_nx_ny_(const unsigned int *nnx,
                                const unsigned int *nny);

      /**
       * Allocate memory, must be called after set nx/ny.
       */
      void fastscape_setup_();

      /**
       * Set FastScape boundary conditions.
       */
      void fastscape_set_bc_(const unsigned int *jbc);

      /**
       * Set FastScape timestep. This will vary based on the ASPECT timestep.
       */
      void fastscape_set_dt_(const double *dtt);

      /**
       * Initialize FastScape topography.
       */
      void fastscape_init_h_(double *hp);

      /**
       * Initialize FastScape silt fraction during a restart.
       */
      void fastscape_init_f_(double *sf);

      /**
       * Set FastScape erosional parameters on land. These parameters will apply to the stream power law (SPL)
       * and hillslope diffusion for basement and sediment. This can be set between timesteps.
       */
      void fastscape_set_erosional_parameters_(double *kkf,
                                               const double *kkfsed,
                                               const double *mm,
                                               const double *nnn,
                                               double *kkd,
                                               const double *kkdsed,
                                               const double *gg1,
                                               const double *gg2,
                                               const double *pp);

      /**
       * Set FastScape marine erosional parameters. This can be set between timesteps.
       */
      void fastscape_set_marine_parameters_(const double *sl,
                                            const double *p1,
                                            const double *p2,
                                            const double *z1,
                                            const double *z2,
                                            const double *r,
                                            const double *l,
                                            const double *kds1,
                                            const double *kds2);

      /**
       * Set advection velocities for FastScape. This can be set between timesteps.
       */
      void fastscape_set_v_(double *ux,
                            double *uy);

      /**
       * Set FastScape uplift rate. This can be set between timesteps.
       */
      void fastscape_set_u_(double *up);

      /**
       * Set FastScape topography. This can be set between timesteps.
       */
      void fastscape_set_h_(double *hp);

      /**
      * Set FastScape basement. This can be set between timesteps. Sediment within FastScape
      * is considered as the difference between the topography and basement, though this may differ
      * from sediment as seen in ASPECT because the FastScape basement only takes the surface
      * velocities into consideration.
      */
      void fastscape_set_basement_(double *b);

      /**
       * Run FastScape for a single FastScape timestep.
       */
      void fastscape_execute_step_();

      /**
       * Create a .VTK file for the FastScape surface within the FastScape folder of the
       * ASPECT output folder.
       */
      void fastscape_named_vtk_(double *fp,
                                const double *vexp,
                                unsigned int *astep,
                                const char *c,
                                const unsigned int *length);

      /**
       * Copy the current FastScape topography.
       */
      void fastscape_copy_h_(double *hp);

      /**
       * Copy the current FastScape basement.
       */
      void fastscape_copy_basement_(double *b);

      /**
       * Copy the current FastScape silt fraction.
       */
      void fastscape_copy_f_(double *sf);

      /**
       * Copy the current FastScape slopes.
       */
      void fastscape_copy_slope_(double *slopep);

      /**
       * Destroy FastScape.
       */
      void fastscape_destroy_();
    }


    template <int dim>
    FastScape<dim>::~FastScape ()
    {
      // It doesn't seem to matter if this is done on all processors or only on the one that runs
      // FastScape as the destroy function checks if the memory is allocated.
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        fastscape_destroy_();
    }

    template <int dim>
    void
    FastScape<dim>::initialize ()
    {
      CitationInfo::add("fastscape");

      AssertThrow(Plugins::plugin_type_matches<const GeometryModel::Box<dim>>(this->get_geometry_model()),
                  ExcMessage("FastScape can only be run with a box geometry model."));

      const GeometryModel::Box<dim> *geometry
        = dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model());

      // Find the id associated with the top boundary and boundaries that call mesh deformation.
      const types::boundary_id top_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id ("top");
      const std::set<types::boundary_id> mesh_deformation_boundary_ids
        = this->get_mesh_deformation_handler().get_active_mesh_deformation_boundary_indicators();

      // Get the deformation type names called for each boundary.
      std::map<types::boundary_id, std::vector<std::string>> mesh_deformation_boundary_indicators_map
        = this->get_mesh_deformation_handler().get_active_mesh_deformation_names();

      // Loop over each mesh deformation boundary, and make sure FastScape is only called on the surface.
      for (const types::boundary_id id : mesh_deformation_boundary_ids)
        {
          const std::vector<std::string> &names = mesh_deformation_boundary_indicators_map[id];
          for (const auto &name : names)
            {
              if (name == "fastscape")
                AssertThrow(id == top_boundary,
                            ExcMessage("FastScape can only be called on the surface boundary."));
            }
        }

      // Several compositional fields are commonly used in conjunction with the FastScape plugin, i.e.
      // "sediment_age" to track the age of the sediment deposited and "deposition_depth" to track the depth
      // with respect to the unperturbed surface of the model domain. Their values are controlled by setting
      // boundary conditions on the top boundary that is deformed by FastScape. While it is useful to track these
      // fields, they are not needed for any function in the FastScape plugin. If they exist however, we need
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

      // Initialize parameters for restarting FastScape
      restart = this->get_parameters().resume_computation;

      // Since we don't open these until we're on one process, we need to check if the
      // restart files exist beforehand.
      if (restart)
        {
          if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            {
              AssertThrow(Utilities::fexists(this->get_output_directory() + "fastscape_elevation_restart.txt"),
                          ExcMessage("Cannot open topography file to restart FastScape."));
              AssertThrow(Utilities::fexists(this->get_output_directory() + "fastscape_basement_restart.txt"),
                          ExcMessage("Cannot open basement file to restart FastScape."));
              AssertThrow(Utilities::fexists(this->get_output_directory() + "fastscape_silt_fraction_restart.txt"),
                          ExcMessage("Cannot open silt fraction file to restart FastScape."));
            }
        }

      // The first entry represents the minimum coordinates of the model domain, the second the model extent.
      for (unsigned int d=0; d<dim; ++d)
        {
          grid_extent[d].first = geometry->get_origin()[d];
          grid_extent[d].second = geometry->get_extents()[d];
        }

      // Get the x and y repetitions used in the parameter file so
      // the FastScape cell size can be properly set.
      const std::array<unsigned int, dim> repetitions = geometry->get_repetitions();

      // Set number of x points, which is generally 1+(FastScape refinement level)^2.
      // The FastScape refinement level is a combination of the maximum ASPECT refinement level
      // at the surface and any additional refinement we want in FastScape. If
      // repetitions are specified we need to adjust the number of points to match what ASPECT has,
      // which can be determined by multiplying the points by the repetitions before adding 1.
      // Finally, if ghost nodes are used we add two additional points on each side.
      const unsigned int ghost_nodes = 2*use_ghost_nodes;
      const unsigned int fastscape_refinement_level = maximum_surface_refinement_level + additional_refinement_levels;
      const unsigned int fastscape_nodes = Utilities::pow(2,fastscape_refinement_level);
      fastscape_nx = fastscape_nodes * repetitions[0] + ghost_nodes + 1;

      // Size of FastScape cell.
      fastscape_dx = (grid_extent[0].second)/(fastscape_nodes * repetitions[0]);

      // FastScape X extent, which is generally ASPECT's extent unless the ghost nodes are used,
      // in which case 2 cells are added on either side.
      fastscape_x_extent = (grid_extent[0].second) + fastscape_dx * ghost_nodes;

      // Sub intervals are 3 less than points, if including the ghost nodes. Otherwise 1 less.
      table_intervals[0] = fastscape_nodes * repetitions[0];
      table_intervals[dim-1] = 1;

      if (dim == 2)
        {
          fastscape_dy = fastscape_dx;
          fastscape_y_extent = std::round(fastscape_y_extent_2d/fastscape_dy)*fastscape_dy + fastscape_dy * ghost_nodes;
          fastscape_ny = 1+fastscape_y_extent/fastscape_dy;
        }
      else
        {
          fastscape_ny = fastscape_nodes * repetitions[1] + ghost_nodes + 1;
          fastscape_dy = (grid_extent[1].second)/(fastscape_nodes * repetitions[1]);
          table_intervals[1] = fastscape_nodes * repetitions[1];
          fastscape_y_extent = (grid_extent[1].second) + fastscape_dy * ghost_nodes;
        }

      // Create a folder for the FastScape visualization files.
      Utilities::create_directory (this->get_output_directory() + "fastscape/",
                                   this->get_mpi_communicator(),
                                   false);

      last_output_time = 0;
    }


    template <int dim>
    void
    FastScape<dim>::compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                             AffineConstraints<double> &mesh_velocity_constraints,
                                                             const std::set<types::boundary_id> &boundary_ids) const
    {

      // Because there is no increase in time during timestep 0, we return and only
      // initialize and run FastScape from timestep 1 and on.
      if (this->get_timestep_number() == 0)
        return;

      TimerOutput::Scope timer_section(this->get_computing_timer(), "FastScape plugin");

      const unsigned int current_timestep = this->get_timestep_number ();
      const double aspect_timestep_in_years = this->get_timestep() / year_in_seconds;

      // Find a FastScape timestep that is below our maximum timestep.
      unsigned int fastscape_iterations = fastscape_steps_per_aspect_step;
      double fastscape_timestep_in_years = aspect_timestep_in_years/fastscape_iterations;
      while (fastscape_timestep_in_years>maximum_fastscape_timestep)
        {
          fastscape_iterations *= 2;
          fastscape_timestep_in_years *= 0.5;
        }

      // Vector to hold the velocities that represent the change to the surface.
      const unsigned int fastscape_array_size = fastscape_nx*fastscape_ny;
      std::vector<double> mesh_velocity_z(fastscape_array_size);

      // FastScape requires multiple specially defined and ordered variables sent to its functions. To make
      // the transfer of these down to one process easier, we first fill out a vector of local_aspect_values,
      // then when we get down to one process we use these local_aspect_values to fill the double arrays
      // in the order needed for FastScape.
      std::vector<std::vector<double>> local_aspect_values = get_aspect_values();

      // Run FastScape on single process.
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          // Initialize the variables that will be sent to FastScape.
          // Elevation is initialized at a very high number so that we can later check that all points
          // received data from ASPECT, and if not throw an assert.
          std::vector<double> elevation(fastscape_array_size, std::numeric_limits<double>::max());
          std::vector<double> velocity_x(fastscape_array_size);
          std::vector<double> velocity_y(fastscape_array_size);
          std::vector<double> velocity_z(fastscape_array_size);
          std::vector<double> bedrock_river_incision_rate_array(fastscape_array_size);
          std::vector<double> bedrock_transport_coefficient_array(fastscape_array_size);
          std::vector<double> basement(fastscape_array_size);
          std::vector<double> silt_fraction(fastscape_array_size);
          std::vector<double> elevation_old(fastscape_array_size);

          fill_fastscape_arrays(elevation,
                                bedrock_transport_coefficient_array,
                                bedrock_river_incision_rate_array,
                                velocity_x,
                                velocity_y,
                                velocity_z,
                                local_aspect_values);

          if (current_timestep == 1 || restart)
            {
              this->get_pcout() << "   Initializing FastScape... " << (1+maximum_surface_refinement_level+additional_refinement_levels) <<
                                " levels, cell size: " << fastscape_dx << " m." << std::endl;

              // Set ghost nodes before initializing.
              if (use_ghost_nodes && !restart)
                set_ghost_nodes(elevation,
                                velocity_x,
                                velocity_y,
                                velocity_z,
                                fastscape_timestep_in_years,
                                true);

              // If we are restarting from a checkpoint, load h values for FastScape instead of using the ASPECT values.
              if (restart)
                {
                  read_restart_files(elevation,
                                     basement,
                                     silt_fraction);

                  restart = false;
                }

              initialize_fastscape(elevation,
                                   basement,
                                   bedrock_transport_coefficient_array,
                                   bedrock_river_incision_rate_array,
                                   silt_fraction);
            }
          else
            {
              // If it isn't the first timestep we ignore initialization and instead copy all height values from FastScape.
              // Generally, we overwrite the topography data from ASPECT as FastScape may be at a higher resolution. However,
              // if we are not using FastScape to advect then we do not want to do this and instead use the ASPECT values.
              if (fastscape_advection_uplift)
                fastscape_copy_h_(elevation.data());
            }

          // Find the appropriate sediment rain based off the time interval.
          const double time_in_years = this->get_time() / year_in_seconds;
          auto it = std::lower_bound(sediment_rain_times.begin(), sediment_rain_times.end(), time_in_years);
          const unsigned int inds = std::distance(sediment_rain_times.begin(), it);
          const double sediment_rain = sediment_rain_rates[inds];

          // Keep initial h values so we can calculate velocity later.
          // In the first timestep, h will be given from other processes.
          // In later timesteps, we copy h directly from FastScape.
          std::mt19937 random_number_generator(fastscape_seed);
          std::uniform_real_distribution<double> random_distribution(-noise_elevation,noise_elevation);
          for (unsigned int i=0; i<fastscape_array_size; ++i)
            {
              elevation_old[i] = elevation[i];

              // Initialize random noise after elevation_old is set, so ASPECT sees this initial topography change.
              // Changing boundary height directly on a fixed FastScape boundary causes reproducibility issues,
              // as such we do not add noise to the boundaries regardless of whether they are ghost nodes
              // or not. However, the boundaries can be changed using the uplift velocity and not cause
              // these issues.
              // TODO: Should this be done through velocities instead of a flat height change?
              if (!is_ghost_node(i,true))
                {
                  if (current_timestep == 1)
                    {
                      // + or - topography based on the initial noise magnitude.
                      const double elevation_seed = random_distribution(random_number_generator);
                      elevation[i] = elevation[i] + elevation_seed;
                    }

                  // Here we add the sediment rain (m/yr) as a flat increase in height.
                  // This is done because adding it as an uplift rate would affect the basement.
                  if (sediment_rain > 0 && use_marine_component)
                    {
                      // Only apply sediment rain to areas below sea level.
                      if (elevation[i] < sea_level)
                        {
                          // If the rain would put us above sea level, set height to sea level.
                          if (elevation[i] + sediment_rain*aspect_timestep_in_years > sea_level)
                            elevation[i] = sea_level;
                          else
                            elevation[i] = std::min(sea_level,elevation[i] + sediment_rain*aspect_timestep_in_years);
                        }
                    }
                }
            }

          // The ghost nodes are added as a single layer of points surrounding the entire model.
          // For example, if ASPECT's surface mesh is a 2D surface that is 3x3 (nx x ny) points,
          // FastScape will be set as a 2D 5x5 point surface. On return to ASPECT, the outer ghost nodes
          // will be ignored, and ASPECT will see only the inner 3x3 surface of FastScape.
          if (use_ghost_nodes)
            set_ghost_nodes(elevation,
                            velocity_x,
                            velocity_y,
                            velocity_z,
                            fastscape_timestep_in_years,
                            false);

          // If specified, apply the orographic controls to the FastScape model.
          if (use_orographic_controls)
            apply_orographic_controls(elevation,
                                      bedrock_transport_coefficient_array,
                                      bedrock_river_incision_rate_array);

          // Set velocity components.
          if (fastscape_advection_uplift)
            {
              fastscape_set_u_(velocity_z.data());
              fastscape_set_v_(velocity_x.data(),
                               velocity_y.data());
            }

          // Set h to new values, and erosional parameters if there have been changes.
          fastscape_set_h_(elevation.data());

          fastscape_set_erosional_parameters_(bedrock_river_incision_rate_array.data(),
                                              &sediment_river_incision_rate,
                                              &drainage_area_exponent_m,
                                              &slope_exponent_n,
                                              bedrock_transport_coefficient_array.data(),
                                              &sediment_transport_coefficient,
                                              &bedrock_deposition_g,
                                              &sediment_deposition_g,
                                              &slope_exponent_p);

          // Find timestep size, run FastScape, and make visualizations.
          execute_fastscape(elevation,
                            bedrock_transport_coefficient_array,
                            velocity_x,
                            velocity_y,
                            velocity_z,
                            fastscape_timestep_in_years,
                            fastscape_iterations);

          // Write a file to store h, b & step for restarting.
          // TODO: It would be good to roll this into the general ASPECT checkpointing,
          // and when we do this needs to be changed.
          if (((this->get_parameters().checkpoint_time_secs == 0) &&
               (this->get_parameters().checkpoint_steps > 0) &&
               ((current_timestep + 1) % this->get_parameters().checkpoint_steps == 0)) ||
              (this->get_time() == this->get_end_time() && this->get_timestepping_manager().need_checkpoint_on_terminate()))
            {
              save_restart_files(elevation,
                                 basement,
                                 silt_fraction);
            }

          // Find out our velocities from the change in height.
          // Where mesh_velocity_z is a vector of array size that exists on all processes.
          for (unsigned int i=0; i<fastscape_array_size; ++i)
            {
              mesh_velocity_z[i] = (elevation[i] - elevation_old[i])/aspect_timestep_in_years;
            }

          Utilities::MPI::broadcast(this->get_mpi_communicator(), mesh_velocity_z, 0);
        }
      else
        {
          for (unsigned int i=0; i<local_aspect_values.size(); ++i)
            MPI_Ssend(&local_aspect_values[i][0], local_aspect_values[1].size(), MPI_DOUBLE, 0, 42, this->get_mpi_communicator());

          // Check whether the FastScape mesh was filled with data.
          const bool fastscape_mesh_filled = Utilities::MPI::broadcast (this->get_mpi_communicator(), true, 0);
          if (fastscape_mesh_filled != true)
            throw aspect::QuietException();

          // This is called solely so we can set the timer and will return immediately.
          execute_fastscape(mesh_velocity_z,
                            mesh_velocity_z,
                            mesh_velocity_z,
                            mesh_velocity_z,
                            mesh_velocity_z,
                            aspect_timestep_in_years,
                            fastscape_steps_per_aspect_step);

          mesh_velocity_z = Utilities::MPI::broadcast(this->get_mpi_communicator(), mesh_velocity_z, 0);
        }

      // Get the sizes needed for a data table of the mesh velocities.
      TableIndices<dim> size_idx;
      for (unsigned int d=0; d<dim; ++d)
        {
          size_idx[d] = table_intervals[d]+1;
        }

      // Initialize a table to hold all velocity values that will be interpolated back to ASPECT.
      const Table<dim,double> velocity_table = fill_data_table(mesh_velocity_z, size_idx, fastscape_nx, fastscape_ny);

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
    FastScape<dim>::get_aspect_values() const
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
                    // Subtract the origin point so that it corresponds to an origin of 0,0 in FastScape.
                    const double indx = 1+use_ghost_nodes+(vertex(0) - grid_extent[0].first)/fastscape_dx;

                    // The quadrature rule is created so that there are enough interpolation points in the
                    // lowest resolved ASPECT surface cell to fill out the FastScape mesh. However, as the
                    // same rule is used for all cell sizes, higher resolution areas will have interpolation
                    // points that do not correspond to a FastScape node. In which case, indx will not be a
                    // whole number and we can ignore the point.
                    if (std::abs(indx - std::round(indx)) >= node_tolerance)
                      continue;


                    // If we're in 2D, we want to take the values and apply them to every row of X points.
                    if (dim == 2)
                      {
                        for (unsigned int ys=0; ys<fastscape_ny; ++ys)
                          {
                            // FastScape indexes from 1 to n, starting at X and Y = 0, and increases
                            // across the X row. At the end of the row, it jumps back to X = 0
                            // and up to the next X row in increasing Y direction. We track
                            // this to correctly place the variables later on.
                            // Nx*ys effectively tells us what row we are in
                            // and then indx tells us what position in that row.
                            const double index = std::round(indx)+fastscape_nx*ys;

                            local_aspect_values[0].push_back(vertex(dim-1) - grid_extent[dim-1].second);
                            local_aspect_values[1].push_back(index-1);

                            for (unsigned int d=0; d<dim; ++d)
                              {
                                // Always convert to m/yr for FastScape
                                local_aspect_values[2+d].push_back(vel[corner][d]*year_in_seconds);
                              }
                          }
                      }
                    // 3D case
                    else
                      {
                        // Because indy only gives us the row we're in, we don't need to add 2 for the ghost node.
                        const double indy = 1+use_ghost_nodes+(vertex(1) - grid_extent[1].first)/fastscape_dy;

                        if (std::abs(indy - std::round(indy)) >= node_tolerance)
                          continue;

                        const double index = std::round((indy-1))*fastscape_nx+std::round(indx);

                        local_aspect_values[0].push_back(vertex(dim-1) - grid_extent[dim-1].second);   //z component
                        local_aspect_values[1].push_back(index-1);

                        for (unsigned int d=0; d<dim; ++d)
                          {
                            local_aspect_values[2+d].push_back(vel[corner][d]*year_in_seconds);
                          }
                      }
                  }
              }

      return local_aspect_values;
    }


    template <int dim>
    void FastScape<dim>::fill_fastscape_arrays(std::vector<double> &elevation,
                                               std::vector<double> &bedrock_transport_coefficient_array,
                                               std::vector<double> &bedrock_river_incision_rate_array,
                                               std::vector<double> &velocity_x,
                                               std::vector<double> &velocity_y,
                                               std::vector<double> &velocity_z,
                                               std::vector<std::vector<double>> &local_aspect_values) const
    {
      for (unsigned int i=0; i<local_aspect_values[1].size(); ++i)
        {

          unsigned int index = local_aspect_values[1][i];
          elevation[index] = local_aspect_values[0][i];
          velocity_x[index] = local_aspect_values[2][i];
          velocity_z[index] = local_aspect_values[dim+1][i];

          if (dim == 2)
            velocity_y[index] = 0;
          else
            velocity_y[index] = local_aspect_values[3][i];
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
              const unsigned int index = local_aspect_values[1][i];
              elevation[index] = local_aspect_values[0][i];
              velocity_x[index] = local_aspect_values[2][i];
              velocity_z[index] = local_aspect_values[dim+1][i];

              // In 2D there are no y velocities, so we set them to zero.
              if (dim == 2 )
                velocity_y[index] = 0;
              else
                velocity_y[index] = local_aspect_values[3][i];
            }
        }

      // Initialize the bedrock river incision rate and transport coefficient,
      // and check that there are no empty mesh points due to
      // an improperly set maximum_surface_refinement_level, additional_refinement_levels,
      // and surface_refinement_difference
      bool fastscape_mesh_filled = true;
      const unsigned int fastscape_array_size = fastscape_nx*fastscape_ny;
      for (unsigned int i=0; i<fastscape_array_size; ++i)
        {
          bedrock_river_incision_rate_array[i] = bedrock_river_incision_rate;
          bedrock_transport_coefficient_array[i] = bedrock_transport_coefficient;

          // If this is a boundary node that is a ghost node then ignore that it
          // has not filled yet as the ghost nodes haven't been set.
          if (elevation[i] == std::numeric_limits<double>::max() && !is_ghost_node(i,false))
            fastscape_mesh_filled = false;
        }

      Utilities::MPI::broadcast(this->get_mpi_communicator(), fastscape_mesh_filled, 0);
      AssertThrow (fastscape_mesh_filled == true,
                   ExcMessage("The FastScape mesh is missing data. A likely cause for this is that the "
                              "maximum surface refinement or surface refinement difference are improperly set."));
    }


    template <int dim>
    void FastScape<dim>::initialize_fastscape(std::vector<double> &elevation,
                                              std::vector<double> &basement,
                                              std::vector<double> &bedrock_transport_coefficient_array,
                                              std::vector<double> &bedrock_river_incision_rate_array,
                                              std::vector<double> &silt_fraction) const
    {
      Assert (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0, ExcInternalError());

      const unsigned int current_timestep = this->get_timestep_number ();

      // Initialize FastScape with grid and extent.
      fastscape_init_();
      fastscape_set_nx_ny_(&fastscape_nx,
                           &fastscape_ny);
      fastscape_setup_();
      fastscape_set_xl_yl_(&fastscape_x_extent,
                           &fastscape_y_extent);

      // Set boundary conditions
      fastscape_set_bc_(&fastscape_boundary_conditions);

      // Initialize topography
      fastscape_init_h_(elevation.data());

      // Set erosional parameters.
      fastscape_set_erosional_parameters_(bedrock_river_incision_rate_array.data(),
                                          &sediment_river_incision_rate,
                                          &drainage_area_exponent_m,
                                          &slope_exponent_n,
                                          bedrock_transport_coefficient_array.data(),
                                          &sediment_transport_coefficient,
                                          &bedrock_deposition_g,
                                          &sediment_deposition_g,
                                          &slope_exponent_p);

      if (use_marine_component)
        fastscape_set_marine_parameters_(&sea_level,
                                         &sand_surface_porosity,
                                         &silt_surface_porosity,
                                         &sand_efold_depth,
                                         &silt_efold_depth,
                                         &sand_silt_ratio,
                                         &sand_silt_averaging_depth,
                                         &sand_transport_coefficient,
                                         &silt_transport_coefficient);

      // Only set the basement and silt_fraction if it's a restart
      if (current_timestep != 1)
        {
          fastscape_set_basement_(basement.data());
          if (use_marine_component)
            fastscape_init_f_(silt_fraction.data());
        }

    }


    template <int dim>
    void FastScape<dim>::execute_fastscape(std::vector<double> &elevation,
                                           std::vector<double> &extra_vtk_field,
                                           std::vector<double> &velocity_x,
                                           std::vector<double> &velocity_y,
                                           std::vector<double> &velocity_z,
                                           const double &fastscape_timestep_in_years,
                                           const unsigned int &fastscape_iterations) const
    {
      TimerOutput::Scope timer_section(this->get_computing_timer(), "Execute FastScape");
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) != 0)
        return;

      // Because on the first timestep we will create an initial VTK file before running FastScape
      // and a second after, we first set the visualization step to zero.
      unsigned int visualization_step = 0;
      const unsigned int current_timestep = this->get_timestep_number ();
      std::string dirname = (this->get_output_directory() + "fastscape/");
      const char *dirname_char=dirname.c_str();
      const unsigned int dirname_length = dirname.length();

      // Set time step
      fastscape_set_dt_(&fastscape_timestep_in_years);
      this->get_pcout() << "   Executing FastScape... " << (fastscape_iterations) << " timesteps of " << fastscape_timestep_in_years << " years." << std::endl;
      {
        // If it is the first timestep, write an initial VTK file.
        if (current_timestep == 1)
          {
            this->get_pcout() << "      Writing initial VTK..." << std::endl;
            // FastScape by default visualizes a field called HHHHH,
            // and the parameter this shows will be whatever is given as the first
            // position. At the moment it visualizes the bedrock diffusivity.
            fastscape_named_vtk_(extra_vtk_field.data(),
                                 &vexp,
                                 &visualization_step,
                                 dirname_char,
                                 &dirname_length);
          }

        for (unsigned int fastscape_iteration = 0; fastscape_iteration < fastscape_iterations; ++fastscape_iteration)
          {
            fastscape_execute_step_();

            // If we are using the ghost nodes we want to reset them every FastScape timestep.
            if (use_ghost_nodes)
              {
                fastscape_copy_h_(elevation.data());

                set_ghost_nodes(elevation,
                                velocity_x,
                                velocity_y,
                                velocity_z,
                                fastscape_timestep_in_years,
                                false);

                // Set velocity components.
                if (fastscape_advection_uplift)
                  {
                    fastscape_set_u_(velocity_z.data());
                    fastscape_set_v_(velocity_x.data(),
                                     velocity_y.data());
                  }

                // Set h to new values, and erosional parameters if there have been changes.
                fastscape_set_h_(elevation.data());
              }
          }

        // Copy h values.
        fastscape_copy_h_(elevation.data());


        // Determine whether to create a VTK file this timestep.
        bool write_vtk = false;

        if (this->get_time() >= last_output_time + output_interval || this->get_time() == this->get_end_time())
          {
            write_vtk = true;

            if (output_interval > 0)
              {
                // We need to find the last time output was supposed to be written.
                // this is the last_output_time plus the largest positive multiple
                // of output_intervals that passed since then. We need to handle the
                // edge case where last_output_time+output_interval==current_time,
                // we did an output and std::floor sadly rounds to zero. This is done
                // by forcing std::floor to round 1.0-eps to 1.0.
                const double magic = 1.0+2.0*std::numeric_limits<double>::epsilon();
                last_output_time = last_output_time + std::floor((this->get_time()-last_output_time)/output_interval*magic) * output_interval/magic;
              }
          }

        if (write_vtk)
          {
            this->get_pcout() << "      Writing FastScape VTK..." << std::endl;
            visualization_step = current_timestep;
            fastscape_named_vtk_(extra_vtk_field.data(),
                                 &vexp,
                                 &visualization_step,
                                 dirname_char,
                                 &dirname_length);
          }
      }
    }


    template <int dim>
    void FastScape<dim>::apply_orographic_controls(const std::vector<double> &elevation,
                                                   std::vector<double> &bedrock_transport_coefficient_array,
                                                   std::vector<double> &bedrock_river_incision_rate_array) const
    {
      // First for the wind barrier, we find the maximum height and index
      // along each line in the x and y direction.
      // If wind is east or west, we find maximum point for each ny row along x.
      const unsigned int fastscape_array_size = fastscape_nx*fastscape_ny;
      std::vector<std::vector<double>> max_elevation_along_x(2, std::vector<double>(fastscape_ny, 0.0));
      if (wind_direction == 0 || wind_direction == 1)
        {
          for (unsigned int i=0; i<fastscape_ny; ++i)
            {
              for (unsigned int j=0; j<fastscape_nx; ++j)
                {
                  if ( elevation[fastscape_nx*i+j] > max_elevation_along_x[0][i])
                    {
                      // Maximum elevation value along the ny row.
                      max_elevation_along_x[0][i] = elevation[fastscape_nx*i+j];
                      // Location of maximum elevation.
                      max_elevation_along_x[1][i] = j;
                    }
                }
            }
        }

      // If wind is north or south, we find maximum point for each nx row along y.
      std::vector<std::vector<double>> max_elevation_along_y(2, std::vector<double>(fastscape_nx, 0.0));
      if (wind_direction == 2 || wind_direction == 3)
        {
          for (unsigned int i=0; i<fastscape_nx; ++i)
            {
              for (unsigned int j=0; j<fastscape_ny; ++j)
                {
                  if ( elevation[fastscape_nx*j+i] > max_elevation_along_y[0][i])
                    {
                      max_elevation_along_y[0][i] = elevation[fastscape_nx*j+i];
                      max_elevation_along_y[1][i] = j;
                    }
                }
            }
        }

      // Now we loop through all the points again and apply the factors.
      std::vector<double> control_applied(fastscape_array_size, 0);
      for (unsigned int i=0; i<fastscape_ny; ++i)
        {
          // Factor from wind barrier. Apply a switch based off wind direction.
          // Where 0 is wind going to the west, 1 the east, 2 the south, and 3 the north.
          for (unsigned int j=0; j<fastscape_nx; ++j)
            {
              switch (wind_direction)
                {
                  case 0 :
                  {
                    // If we are above the set elevation, and on the correct side based on the wind direction apply
                    // the factor. Apply this regardless of whether or not we stack controls.
                    if ( (max_elevation_along_x[0][i] > wind_barrier_elevation) && (j < max_elevation_along_x[1][i]) )
                      {
                        bedrock_river_incision_rate_array[fastscape_nx*i+j] *= wind_barrier_erosional_factor;
                        bedrock_transport_coefficient_array[fastscape_nx*i+j] *= wind_barrier_erosional_factor;
                        control_applied[fastscape_nx*i+j] = 1;
                      }
                    break;
                  }
                  case 1 :
                  {
                    if ( (max_elevation_along_x[0][i] > wind_barrier_elevation) && (j > max_elevation_along_x[1][i]) )
                      {
                        bedrock_river_incision_rate_array[fastscape_nx*i+j] *= wind_barrier_erosional_factor;
                        bedrock_transport_coefficient_array[fastscape_nx*i+j] *= wind_barrier_erosional_factor;
                        control_applied[fastscape_nx*i+j] = 1;
                      }
                    break;
                  }
                  case 2 :
                  {
                    if ( (max_elevation_along_y[0][j] > wind_barrier_elevation) && (i > max_elevation_along_y[1][j]) )
                      {
                        bedrock_river_incision_rate_array[fastscape_nx*i+j] *= wind_barrier_erosional_factor;
                        bedrock_transport_coefficient_array[fastscape_nx*i+j] *= wind_barrier_erosional_factor;
                        control_applied[fastscape_nx*i+j] = 1;
                      }
                    break;
                  }
                  case 3 :
                  {
                    if ( (max_elevation_along_y[0][j] > wind_barrier_elevation) && (i < max_elevation_along_y[1][j]) )
                      {
                        bedrock_river_incision_rate_array[fastscape_nx*i+j] *= wind_barrier_erosional_factor;
                        bedrock_transport_coefficient_array[fastscape_nx*i+j] *= wind_barrier_erosional_factor;
                        control_applied[fastscape_nx*i+j] = 1;
                      }
                    break;
                  }
                  default :
                    AssertThrow(false, ExcMessage("This does not correspond with a wind direction."));
                    break;
                }

              // If we are above the flat elevation and stack controls, apply the flat elevation factor. If we are not
              // stacking controls, apply the factor if the wind barrier was not applied to this point.
              if (elevation[fastscape_nx*i+j] > flat_elevation)
                {
                  if ( stack_controls==true || !stack_controls && (control_applied[fastscape_nx*i+j]==0) )
                    {
                      bedrock_river_incision_rate_array[fastscape_nx*i+j] *= flat_erosional_factor;
                      bedrock_transport_coefficient_array[fastscape_nx*i+j] *= flat_erosional_factor;
                    }
                  // If we are not stacking controls and the wind barrier was applied to this point, only
                  // switch to this control if the factor is greater.
                  else if ( stack_controls==false && (control_applied[fastscape_nx*i+j]==1) && (flat_erosional_factor > wind_barrier_erosional_factor) )
                    {
                      if ( wind_barrier_erosional_factor != 0)
                        {
                          bedrock_river_incision_rate_array[fastscape_nx*i+j] = (bedrock_river_incision_rate_array[fastscape_nx*i+j]/wind_barrier_erosional_factor)*flat_erosional_factor;
                          bedrock_transport_coefficient_array[fastscape_nx*i+j] = (bedrock_transport_coefficient_array[fastscape_nx*i+j]/wind_barrier_erosional_factor)*flat_erosional_factor;
                        }
                      // If a wind barrier factor of zero was applied for some reason, we set it back to the default
                      // and apply the flat_erosional_factor.
                      else
                        {
                          bedrock_river_incision_rate_array[fastscape_nx*i+j] = bedrock_river_incision_rate*flat_erosional_factor;
                          bedrock_transport_coefficient_array[fastscape_nx*i+j] = bedrock_transport_coefficient*flat_erosional_factor;
                        }
                    }
                }
            }
        }
    }


    template <int dim>
    void FastScape<dim>::set_ghost_nodes(std::vector<double> &elevation,
                                         std::vector<double> &velocity_x,
                                         std::vector<double> &velocity_y,
                                         std::vector<double> &velocity_z,
                                         const double &fastscape_timestep_in_years,
                                         const bool init) const
    {
      // Copy the slopes at each point, this will be used to set an H
      // at the ghost nodes if a boundary mass flux is given.
      const unsigned int fastscape_array_size = fastscape_nx*fastscape_ny;
      std::vector<double> slopep(fastscape_array_size);

      if (!init)
        fastscape_copy_slope_(slopep.data());

      // Here we set the ghost nodes at the left and right boundaries. In most cases,
      // this involves setting the node to the same values of v and h as the inward node.
      // With the inward node being above or below for the bottom and top rows of ghost nodes,
      // or to the left and right for the right and left columns of ghost nodes.
      for (unsigned int j=0; j<fastscape_ny; ++j)
        {
          // Nx*j will give us the row we're in, and one is subtracted as FastScape starts from 1 not zero.
          // If we're on the left, the multiple of the row will always represent the first node.
          // Subtracting one from the row above this gives us the last node of the previous row.
          const unsigned int index_left = fastscape_nx*j;
          const unsigned int index_right = fastscape_nx*(j+1)-1;
          double slope = 0;

          // Here we set the ghost nodes to the value of the nodes next to them, where for the left we
          // add one to go to the node to the right, and for the right side
          // we subtract one to go to the inner node to the left.
          // For velocities, this is done whether it is the call before
          // initialization or not.
          if (left == 0 || !use_fixed_erosional_base)
            {
              velocity_z[index_left] = velocity_z[index_left+1];
              velocity_y[index_left] = velocity_y[index_left+1];
              velocity_x[index_left] = velocity_x[index_left+1];
            }

          if (right == 0 || !use_fixed_erosional_base)
            {
              velocity_z[index_right] = velocity_z[index_right-1];
              velocity_y[index_right] = velocity_y[index_right-1];
              velocity_x[index_right] = velocity_x[index_right-1];
            }

          if (init)
            {
              // We need to set h on initialization. If we use a fixed base level we set it to
              // the user specified base level. If it is open or we do not have a fixed base, we
              // will set it equal to the node next to it, and finally adjust based on user-defined
              // influx if necessary. FastScape calculates the slope by looking at all
              // nodes surrounding the point so we need to consider the slope over 2 dx.
              slope = left_flux/bedrock_transport_coefficient;
              if (left == 1 && use_fixed_erosional_base)
                elevation[index_left] = h_erosional_base;
              else
                elevation[index_left] = elevation[index_left+1] + slope*2*fastscape_dx;

              slope = right_flux/bedrock_transport_coefficient;
              if (right == 1 && use_fixed_erosional_base)
                elevation[index_right] = h_erosional_base;
              else
                elevation[index_right] = elevation[index_right-1] + slope*2*fastscape_dx;
            }

          if (left == 0 && !init)
            {
              // If it is not the initialization step, we set the h value for open boundaries depending on
              // whether or not there is prescribed influx.
              // If we have flux through a boundary, we need to update the height to keep the correct slope.
              // Because the corner nodes always show a slope of zero, this will update them according to
              // the closest non-ghost node. E.g. if we're at a corner node, look instead up a row and inward.
              // If this is no flux, we set the node to the one next to it.
              if (left_flux > 0)
                {
                  slope = 0;
                  if (j == 0)
                    slope = left_flux / bedrock_transport_coefficient - std::tan(slopep[index_left + fastscape_nx + 1] * numbers::PI / 180.);
                  else if (j == (fastscape_ny - 1))
                    slope = left_flux / bedrock_transport_coefficient - std::tan(slopep[index_left - fastscape_nx + 1] * numbers::PI / 180.);
                  else
                    slope = left_flux / bedrock_transport_coefficient - std::tan(slopep[index_left + 1] * numbers::PI / 180.);

                  elevation[index_left] = elevation[index_left] + slope * 2 * fastscape_dx;
                }
              else
                elevation[index_left] = elevation[index_left + 1];
            }

          if (right == 0 && !init)
            {

              if (right_flux > 0)
                {
                  slope = 0;
                  if (j == 0)
                    slope = right_flux / bedrock_transport_coefficient - std::tan(slopep[index_right + fastscape_nx - 1] * numbers::PI / 180.);
                  else if (j == (fastscape_ny - 1))
                    slope = right_flux / bedrock_transport_coefficient - std::tan(slopep[index_right - fastscape_nx - 1] * numbers::PI / 180.);
                  else
                    slope = right_flux / bedrock_transport_coefficient - std::tan(slopep[index_right - 1] * numbers::PI / 180.);

                  elevation[index_right] = elevation[index_right] + slope * 2 * fastscape_dx;
                }
              else
                elevation[index_right] = elevation[index_right - 1];


            }

          // If the boundaries are periodic, then we look at the velocities on both sides of the
          // model, and set the ghost node according to the direction of flow. As FastScape will
          // receive all velocities it will have a direction, and we only need to look at the (non-ghost)
          // nodes directly to the left and right.
          if (left == 0 && right == 0 || leftright_ghost_nodes_periodic == true)
            {
              // First we assume that flow is going to the left.
              unsigned int side = index_left;
              unsigned int op_side = index_right;

              // Indexing depending on which side the ghost node is being set to.
              int jj = 1;

              // If nodes on both sides are going the same direction, then set the respective
              // ghost nodes to equal these sides. By doing this, the ghost nodes at the opposite
              // side of flow will work as a mirror mimicking what is happening on the other side.
              if (velocity_x[index_right-1] > 0 && velocity_x[index_left+1] >= 0)
                {
                  side = index_right;
                  op_side = index_left;
                  jj = -1;
                }
              else if (velocity_x[index_right-1] <= 0 && velocity_x[index_left+1] < 0)
                {
                  side = index_left;
                  op_side = index_right;
                  jj = 1;
                }
              else
                continue;

              // Now set the nodes for periodic boundaries. As an example, assume we have 9 FastScape nodes in x:
              //
              // 0 - 1 - 2 - 3 - 4 - 5 - 6 - 7 - 8
              //
              // Of these 9 nodes, 0 and 8 are ghost nodes and 1 and 7 are the periodic ASPECT boundaries.
              // If we assume that the horizontal ASPECT direction of travel is towards node 7, then we would
              // set the ghost node 8 velocities and heights to that of node 2, ghost node 0 to node 6, and
              // ASPECT boundary node 1 to ASPECT boundary node 7. E.g., based on the FastScape values for
              // vx, vy, vz, and elevation, the nodes could be rewritten as:
              //
              // 6 - 7 - 2 - 3 - 4 - 5 - 6 - 7 - 2
              //
              // This makes it so that effectively both periodic ASPECT boundaries see the same
              // topography on either side of them to try and make sure they experience the same
              // amount of diffusion and SPL.
              velocity_x[index_right] = velocity_x[index_left+2];
              velocity_y[index_right] = velocity_y[index_left+2];
              velocity_z[index_right] = velocity_z[index_left+2] + (elevation[index_left+2] - elevation[index_right])/fastscape_timestep_in_years;

              velocity_x[index_left] = velocity_x[index_right-2];
              velocity_y[index_left] = velocity_y[index_right-2];
              velocity_z[index_left] = velocity_z[index_right-2] + (elevation[index_right-2] - elevation[index_left])/fastscape_timestep_in_years;

              // Set opposing ASPECT boundary so it's periodic.
              elevation[op_side-jj] = elevation[side+jj];
              velocity_x[op_side-jj] = velocity_x[side+jj];
              velocity_y[op_side-jj] = velocity_y[side+jj];
              velocity_z[op_side-jj] = velocity_z[side+jj];

            }
        }

      // Now do the same for the top and bottom ghost nodes.
      for (unsigned int j=0; j<fastscape_nx; ++j)
        {
          // The bottom row indexes are 0 to nx-1.
          const unsigned int index_bot = j;

          // Nx multiplied by (total rows - 1) gives us the start of
          // the top row, and j gives the position in the row.
          const unsigned int index_top = fastscape_nx*(fastscape_ny-1)+j;
          double slope = 0;

          if (top == 0 || !use_fixed_erosional_base)
            {
              velocity_z[index_top] = velocity_z[index_top-fastscape_nx];
              velocity_y[index_top] = velocity_y[index_top-fastscape_nx];
              velocity_x[index_top] = velocity_x[index_top-fastscape_nx];
            }

          if (bottom ==0 || !use_fixed_erosional_base)
            {
              velocity_z[index_bot] = velocity_z[index_bot+fastscape_nx];
              velocity_y[index_bot] = velocity_y[index_bot+fastscape_nx];
              velocity_x[index_bot] = velocity_x[index_bot+fastscape_nx];
            }

          if (init)
            {
              slope = top_flux / bedrock_transport_coefficient;
              if (top == 1 && use_fixed_erosional_base)
                elevation[index_top] = h_erosional_base;
              else
                elevation[index_top] = elevation[index_top-fastscape_nx] + slope*2*fastscape_dx;

              slope = bottom_flux / bedrock_transport_coefficient;
              if (bottom == 1 && use_fixed_erosional_base)
                elevation[index_bot] = h_erosional_base;
              else
                elevation[index_bot] = elevation[index_bot + fastscape_nx] + slope*2*fastscape_dx;
            }

          if (top == 0 && !init)
            {
              if (top_flux > 0)
                {
                  slope = 0;
                  if (j == 0)
                    slope = top_flux / bedrock_transport_coefficient - std::tan(slopep[index_top - fastscape_nx + 1] * numbers::PI / 180.);
                  else if (j == (fastscape_nx - 1))
                    slope = top_flux / bedrock_transport_coefficient - std::tan(slopep[index_top - fastscape_nx - 1] * numbers::PI / 180.);
                  else
                    slope = top_flux / bedrock_transport_coefficient - std::tan(slopep[index_top - fastscape_nx] * numbers::PI / 180.);

                  elevation[index_top] = elevation[index_top] + slope * 2 * fastscape_dx;
                }
              else
                elevation[index_top] = elevation[index_top - fastscape_nx];
            }

          if (bottom == 0 && !init)
            {
              if (left_flux > 0)
                {
                  slope = 0;
                  if (j == 0)
                    slope = bottom_flux / bedrock_transport_coefficient - std::tan(slopep[index_bot + fastscape_nx + 1] * numbers::PI / 180.);
                  else if (j == (fastscape_nx - 1))
                    slope = bottom_flux / bedrock_transport_coefficient - std::tan(slopep[index_bot + fastscape_nx - 1] * numbers::PI / 180.);
                  else
                    slope = bottom_flux / bedrock_transport_coefficient - std::tan(slopep[index_bot + fastscape_nx] * numbers::PI / 180.);

                  elevation[index_bot] = elevation[index_bot] + slope * 2 * fastscape_dx;
                }
              else
                elevation[index_bot] = elevation[index_bot + fastscape_nx];
            }

          if (bottom == 0 && top == 0 || topbottom_ghost_nodes_periodic == true)
            {
              unsigned int side = index_bot;
              unsigned int op_side = index_top;
              int jj = fastscape_nx;

              if (velocity_y[index_bot+fastscape_nx-1] > 0 && velocity_y[index_top-fastscape_nx-1] >= 0)
                {
                  side = index_top;
                  op_side = index_bot;
                  jj = -fastscape_nx;
                }
              else if (velocity_y[index_bot+fastscape_nx-1] <= 0 && velocity_y[index_top-fastscape_nx-1] < 0)
                {
                  side = index_bot;
                  op_side = index_top;
                  jj = fastscape_nx;
                }
              else
                continue;

              // Set top ghost node
              velocity_x[index_top] = velocity_x[index_bot + 2*fastscape_nx];
              velocity_y[index_top] = velocity_y[index_bot + 2*fastscape_nx];
              velocity_z[index_top] = velocity_z[index_bot + 2*fastscape_nx] + (elevation[index_bot + 2*fastscape_nx] - elevation[index_top])/fastscape_timestep_in_years;

              // Set bottom ghost node
              velocity_x[index_bot] = velocity_x[index_top - 2*fastscape_nx];
              velocity_y[index_bot] = velocity_y[index_top - 2*fastscape_nx];
              velocity_z[index_bot] = velocity_z[index_top - 2*fastscape_nx] + (elevation[index_top - 2*fastscape_nx] - elevation[index_bot])/fastscape_timestep_in_years;

              // Set opposing ASPECT boundary so it's periodic.
              elevation[op_side-jj] = elevation[side+jj];
              velocity_x[op_side-jj] = velocity_x[side+jj];
              velocity_y[op_side-jj] = velocity_y[side+jj];
              velocity_z[op_side-jj] = velocity_z[side+jj];
            }
        }
    }

    template <int dim>
    bool FastScape<dim>::is_ghost_node(const unsigned int &index,
                                       const bool &exclude_boundaries) const
    {
      if (use_ghost_nodes == false && exclude_boundaries == false)
        return false;

      const unsigned int row = index / fastscape_nx; // Calculate the row index
      const unsigned int col = index % fastscape_nx; // Calculate the column index

      // If we are at a boundary node and ghost nodes are enabled
      // or we are excluding the boundaries then return true.
      if (row == 0 || row == fastscape_ny-1 || col == 0 || col == fastscape_nx-1)
        return true;
      else
        return false;
    }


    template <int dim>
    Table<dim,double>
    FastScape<dim>::fill_data_table(std::vector<double> &values,
                                    TableIndices<dim> &size_idx,
                                    const unsigned int &fastscape_nx,
                                    const unsigned int &fastscape_ny) const
    {
      // Create data table based off of the given size.
      Table<dim,double> data_table;
      data_table.TableBase<dim,double>::reinit(size_idx);
      TableIndices<dim> idx;

      // Loop through the data table and fill it with the velocities from FastScape.
      if (dim == 2)
        {
          std::vector<double> values_2d(fastscape_nx);

          for (unsigned int x=use_ghost_nodes; x<(fastscape_nx-use_ghost_nodes); ++x)
            {
              // If we do not average the values, then use a slice near the center.
              if (!average_out_of_plane_surface_topography)
                {
                  const unsigned int index = x+fastscape_nx*(std::round((fastscape_ny-use_ghost_nodes)/2));

                  // If we are using the ghost nodes, then the x value locations need to be shifted back 1
                  // e.g., given a 4x4 mesh an index of 5 would correspond to an x of 1 and y of 1 in the loop,
                  // but should correspond to 0,0 for ASPECT.
                  values_2d[x-use_ghost_nodes] = values[index];
                }
              // Here we use average velocities across the y nodes, excluding the ghost nodes (top and bottom row).
              // Note: If ghost nodes are turned off, boundary effects may influence this.
              else
                {
                  for (unsigned int y=use_ghost_nodes; y<(fastscape_ny-use_ghost_nodes); ++y)
                    {
                      const unsigned int index = x+fastscape_nx*y;
                      values_2d[x-use_ghost_nodes] += values[index];
                    }
                  values_2d[x-use_ghost_nodes] = values_2d[x-use_ghost_nodes]/(fastscape_ny-2*use_ghost_nodes);
                }
            }

          for (unsigned int x=0; x<data_table.size()[0]; ++x)
            {
              idx[0] = x;

              for (unsigned int y=0; y<(data_table.size()[1]); ++y)
                {
                  idx[1] = y;

                  // Convert back to m/s.
                  data_table(idx) = values_2d[x] / year_in_seconds;
                }
            }
        }
      else
        {
          // Indexes through x, y, and z.
          for (unsigned int x=0; x<data_table.size()[0]; ++x)
            {
              idx[0] = x;

              for (unsigned int y=0; y<data_table.size()[1]; ++y)
                {
                  idx[1] = y;

                  for (unsigned int z=0; z<data_table.size()[2]; ++z)
                    {
                      idx[2] = z;

                      // Convert back to m/s.
                      data_table(idx) = values[(fastscape_nx+1)*use_ghost_nodes+fastscape_nx*y+x] / year_in_seconds;

                    }
                }
            }
        }

      return data_table;
    }



    template <int dim>
    void FastScape<dim>::read_restart_files(std::vector<double> &elevation,
                                            std::vector<double> &basement,
                                            std::vector<double> &silt_fraction) const
    {
      this->get_pcout() << "   Loading FastScape restart file... " << std::endl;

      // Create variables for output directory and restart file
      const unsigned int fastscape_array_size = fastscape_nx*fastscape_ny;
      std::string dirname = this->get_output_directory();
      const std::string restart_filename_elevation = dirname + "fastscape_elevation_restart.txt";
      const std::string restart_filename_basement = dirname + "fastscape_basement_restart.txt";
      const std::string restart_filename_silt_fraction = dirname + "fastscape_silt_fraction_restart.txt";
      const std::string restart_filename_time = dirname + "fastscape_last_output_time.txt";

      // Load in h values.
      std::ifstream in_elevation(restart_filename_elevation);
      AssertThrow (in_elevation, ExcIO());
      {
        unsigned int line = 0;
        while (line < fastscape_array_size)
          {
            in_elevation >> elevation[line];
            line++;
          }
      }

      // Load in b values.
      std::ifstream in_basement(restart_filename_basement);
      AssertThrow (in_basement, ExcIO());
      {
        unsigned int line = 0;
        while (line < fastscape_array_size)
          {
            in_basement >> basement[line];
            line++;
          }
      }

      // Load in silt_fraction values if
      // marine sediment transport and deposition is active.
      if (use_marine_component)
        {
          std::ifstream in_silt_fraction(restart_filename_silt_fraction);
          AssertThrow (in_silt_fraction, ExcIO());

          if (sand_surface_porosity > 0. || silt_surface_porosity > 0.)
            this->get_pcout() << "   Restarting runs with nonzero porosity can lead to a different system after restart. " << std::endl;
          unsigned int line = 0;
          while (line < fastscape_array_size)
            {
              in_silt_fraction >> silt_fraction[line];
              line++;
            }
        }

      // Now load the last output at time of restart.
      // this allows us to correctly track when to call
      // FastScape to make new VTK files.
      std::ifstream in_last_output_time(restart_filename_time);
      AssertThrow (in_last_output_time, ExcIO());
      {
        in_last_output_time >> last_output_time;
      }
    }

    template <int dim>
    void FastScape<dim>::save_restart_files(const std::vector<double> &elevation,
                                            std::vector<double> &basement,
                                            std::vector<double> &silt_fraction) const
    {
      this->get_pcout() << "      Writing FastScape restart file... " << std::endl;

      // Create variables for output directory and restart file
      const unsigned int fastscape_array_size = fastscape_nx*fastscape_ny;
      std::string dirname = this->get_output_directory();
      const std::string restart_filename_elevation = dirname + "fastscape_elevation_restart.txt";
      const std::string restart_filename_basement = dirname + "fastscape_basement_restart.txt";
      const std::string restart_filename_silt_fraction = dirname + "fastscape_silt_fraction_restart.txt";
      const std::string restart_filename_time = dirname + "fastscape_last_output_time.txt";

      std::ofstream out_elevation(restart_filename_elevation);
      std::ofstream out_basement(restart_filename_basement);
      std::ofstream out_silt_fraction(restart_filename_silt_fraction);
      std::ofstream out_last_output_time(restart_filename_time);
      std::stringstream buffer_basement;
      std::stringstream buffer_elevation;
      std::stringstream buffer_silt_fraction;
      std::stringstream buffer_time;

      fastscape_copy_basement_(basement.data());

      // If marine sediment transport and deposition is active,
      // we also need to store the silt fraction.
      if (use_marine_component)
        fastscape_copy_f_(silt_fraction.data());

      out_last_output_time << last_output_time << "\n";

      for (unsigned int i = 0; i < fastscape_array_size; ++i)
        {
          buffer_elevation << elevation[i] << "\n";
          buffer_basement << basement[i] << "\n";
          if (use_marine_component)
            buffer_silt_fraction << silt_fraction[i] << "\n";
        }

      out_elevation << buffer_elevation.str();
      out_basement << buffer_basement.str();
      if (use_marine_component)
        out_silt_fraction << buffer_silt_fraction.str();
    }



    template <int dim>
    bool
    FastScape<dim>::
    needs_surface_stabilization () const
    {
      return true;
    }



    template <int dim>
    void FastScape<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection ("Fastscape");
        {
          prm.declare_entry("Number of fastscape timesteps per aspect timestep", "5",
                            Patterns::Integer(),
                            "Initial number of fastscape time steps per ASPECT timestep, this value will double if"
                            "the FastScape timestep is above the maximum FastScape timestep.");
          prm.declare_entry("Maximum timestep length", "10e3",
                            Patterns::Double(0),
                            "Maximum timestep for FastScape. Units: $\\{yrs}$");
          prm.declare_entry("Vertical exaggeration", "-1",
                            Patterns::Double(),
                            "Vertical exaggeration for FastScape's VTK file. -1 outputs topography, basement, and sealevel.");
          prm.declare_entry("Additional fastscape refinement", "0",
                            Patterns::Integer(),
                            "How many levels above ASPECT FastScape should be refined.");
          prm.declare_entry ("Average out of plane surface topography in 2d", "true",
                             Patterns::Bool (),
                             "If this is set to false, then a 2D model will only consider the "
                             "center slice FastScape gives. If set to true, then ASPECT will"
                             "average the mesh along Y excluding the ghost nodes.");
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
          prm.declare_entry ("Use marine component", "false",
                             Patterns::Bool (),
                             "Flag to use the marine component of FastScape.");
          prm.declare_entry("Y extent in 2d", "100000",
                            Patterns::Double(),
                            "FastScape Y extent when using a 2D ASPECT model. Units: $\\{m}$");
          prm.declare_entry ("Use ghost nodes", "true",
                             Patterns::Bool (),
                             "Flag to use ghost nodes");
          prm.declare_entry ("Uplift and advect with fastscape", "true",
                             Patterns::Bool (),
                             "Flag to use FastScape advection and uplift.");
          prm.declare_entry("Node tolerance", "0.001",
                            Patterns::Double(),
                            "Node tolerance for how close an ASPECT node must be to the FastScape node for the value to be transferred.");
          prm.declare_entry ("Sediment rain rates", "0,0",
                             Patterns::List (Patterns::Double(0)),
                             "Sediment rain rates given as a list 1 greater than the number of sediment rain time intervals. E.g, "
                             " If the time interval is given at 5 Myr, there will be one value for 0-5 Myr model time and a second value "
                             " for 5+ Myr. Units: $\\{m/yr}$");
          prm.declare_entry ("Sediment rain time intervals", "0",
                             Patterns::List (Patterns::Double(0)),
                             "A list of times to change the sediment rain rate. Units: $\\{yrs}$");
          prm.declare_entry("Initial noise magnitude", "5",
                            Patterns::Double(),
                            "Maximum topography change from the initial noise. Units: $\\{m}$");

          prm.enter_subsection ("Boundary conditions");
          {
            prm.declare_entry ("Front", "1",
                               Patterns::Integer (0, 1),
                               "Front (bottom) boundary condition, where 1 is fixed and 0 is reflective.");
            prm.declare_entry ("Right", "1",
                               Patterns::Integer (0, 1),
                               "Right boundary condition, where 1 is fixed and 0 is reflective.");
            prm.declare_entry ("Back", "1",
                               Patterns::Integer (0, 1),
                               "Back (top) boundary condition, where 1 is fixed and 0 is reflective.");
            prm.declare_entry ("Left", "1",
                               Patterns::Integer (0, 1),
                               "Left boundary condition, where 1 is fixed and 0 is reflective.");
            prm.declare_entry("Left mass flux", "0",
                              Patterns::Double(),
                              "Flux per unit length through left boundary. Units: $\\{m^2/yr}$ ");
            prm.declare_entry("Right mass flux", "0",
                              Patterns::Double(),
                              "Flux per unit length through right boundary. Units: $\\{m^2/yr}$ ");
            prm.declare_entry("Back mass flux", "0",
                              Patterns::Double(),
                              "Flux per unit length through back boundary. Units: $\\{m^2/yr}$ ");
            prm.declare_entry("Front mass flux", "0",
                              Patterns::Double(),
                              "Flux per unit length through front boundary. Units: $\\{m^2/yr}$ ");
            prm.declare_entry ("Back front ghost nodes periodic", "false",
                               Patterns::Bool (),
                               "Whether to set the ghost nodes at the FastScape back and front boundary "
                               "to periodic even if 'Back' and 'Front' are set to fixed boundary.");
            prm.declare_entry ("Left right ghost nodes periodic", "false",
                               Patterns::Bool (),
                               "Whether to set the ghost nodes at the FastScape left and right boundary "
                               "to periodic even if 'Left' and 'Right' are set to fixed boundary.");
          }
          prm.leave_subsection();

          prm.enter_subsection ("Erosional parameters");
          {
            prm.declare_entry("Drainage area exponent", "0.4",
                              Patterns::Double(),
                              "Exponent for drainage area.");
            prm.declare_entry("Slope exponent", "1",
                              Patterns::Double(),
                              "The  slope  exponent  for  SPL (n).  Generally  m/n  should  equal  approximately 0.4");
            prm.declare_entry("Multi-direction slope exponent", "1",
                              Patterns::Double(),
                              "Exponent to determine the distribution from the SPL to neighbor nodes, with"
                              "10 being steepest decent and 1 being more varied.");
            prm.declare_entry("Bedrock deposition coefficient", "1",
                              Patterns::Double(),
                              "Deposition coefficient for bedrock.");
            prm.declare_entry("Sediment deposition coefficient", "-1",
                              Patterns::Double(),
                              "Deposition coefficient for sediment, -1 sets this to the same as the bedrock deposition coefficient.");
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
            prm.declare_entry("Orographic elevation control", "2000",
                              Patterns::Integer(),
                              "Above this height, the elevation factor is applied. Units: $\\{m}$");
            prm.declare_entry("Orographic wind barrier height", "500",
                              Patterns::Integer(),
                              "When terrain reaches this height the wind barrier factor is applied. Units: $\\{m}$");
            prm.declare_entry("Elevation factor", "1",
                              Patterns::Double(),
                              "Amount to multiply the bedrock river incision rate nad transport coefficient by past the given orographic elevation control.");
            prm.declare_entry("Wind barrier factor", "1",
                              Patterns::Double(),
                              "Amount to multiply the bedrock river incision rate nad transport coefficient by past given wind barrier height.");
            prm.declare_entry ("Stack orographic controls", "true",
                               Patterns::Bool (),
                               "Whether or not to apply both controls to a point, or only a maximum of one set as the wind barrier.");
            prm.declare_entry ("Flag to use orographic controls", "false",
                               Patterns::Bool (),
                               "Whether or not to apply orographic controls.");
            prm.declare_entry ("Wind direction", "west",
                               Patterns::Selection("east|west|south|north"),
                               "This parameter assumes a wind direction, deciding which side is reduced from the wind barrier.");
            prm.declare_entry ("Use a fixed erosional base level", "false",
                               Patterns::Bool (),
                               "Whether or not to use an erosional base level that differs from sea level. Setting this parameter to "
                               "true will set all ghost nodes of fixed FastScape boundaries to the height you specify in "
                               "'set Erosional base level'. \nThis can make "
                               "sense for a continental model where the model surrounding topography is assumed above sea level, "
                               "e.g. highlands. If the sea level would be used as an erosional base level in this case, all topography "
                               "erodes away with lots of 'sediment volume' lost through the sides of the model. This is mostly "
                               "important, when there are mountains in the middle of the model, while it is less important when there "
                               "is lower relief in the middle of the model. \n"
                               "In the FastScape  visualization files, setting the extra base level may show up as a strong "
                               "slope at the fixed boundaries of the model. However, in the ASPECT visualization files it will not "
                               "show up, as the ghost nodes only exist in FastScape.");
            prm.declare_entry("Erosional base level", "0",
                              Patterns::Double(),
                              "When 'Use a fixed erosional base level' is set to true, all ghost nodes of fixed "
                              "FastScape boundaries where no mass flux is specified by the user (FastScape boundary condition set to 1 "
                              "and 'Left/Right/Bottom/Top mass flux' set to 0) will be fixed to this elevation. The "
                              "reflecting boundaries (FastScape boundary condition set to 0) will not be affected, nor are the "
                              "boundaries where a mass flux is specified. \n"
                              "Units: m");
          }
          prm.leave_subsection();

          prm.enter_subsection ("Marine parameters");
          {
            prm.declare_entry("Sea level", "0",
                              Patterns::Double(),
                              "Sea level relative to the ASPECT surface, where the maximum Z or Y extent in ASPECT is a sea level of zero. Units: $\\{m}$ ");
            prm.declare_entry("Sand porosity", "0.0",
                              Patterns::Double(),
                              "Porosity of sand. ");
            prm.declare_entry("Silt porosity", "0.0",
                              Patterns::Double(),
                              "Porosity of silt. ");
            prm.declare_entry("Sand e-folding depth", "1e3",
                              Patterns::Double(),
                              "E-folding depth for the exponential of the sand porosity law. Units: $\\{m}$");
            prm.declare_entry("Silt e-folding depth", "1e3",
                              Patterns::Double(),
                              "E-folding depth for the exponential of the silt porosity law. Units: $\\{m}$");
            prm.declare_entry("Sand-silt ratio", "0.5",
                              Patterns::Double(),
                              "Ratio of sand to silt for material leaving continent.");
            prm.declare_entry("Depth averaging thickness", "1e2",
                              Patterns::Double(),
                              "Depth averaging for the sand-silt equation. Units: $\\{m}$");
            prm.declare_entry("Sand transport coefficient", "5e2",
                              Patterns::Double(),
                              "Transport coefficient (diffusivity) for sand. Units: $\\{m^2/yr}$");
            prm.declare_entry("Silt transport coefficient", "2.5e2",
                              Patterns::Double(),
                              "Transport coefficient (diffusivity) for silt. Units: $\\{m^2/yr}$ ");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void FastScape<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection("Fastscape");
        {
          fastscape_steps_per_aspect_step = prm.get_integer("Number of fastscape timesteps per aspect timestep");
          maximum_fastscape_timestep = prm.get_double("Maximum timestep length");
          vexp = prm.get_double("Vertical exaggeration");
          additional_refinement_levels = prm.get_integer("Additional fastscape refinement");
          average_out_of_plane_surface_topography = prm.get_bool("Average out of plane surface topography in 2d");
          fastscape_seed = prm.get_integer("Fastscape seed");
          maximum_surface_refinement_level = prm.get_integer("Maximum surface refinement level");
          surface_refinement_difference = prm.get_integer("Surface refinement difference");
          use_marine_component = prm.get_bool("Use marine component");
          fastscape_y_extent_2d = prm.get_double("Y extent in 2d");
          use_ghost_nodes = prm.get_bool("Use ghost nodes");
          fastscape_advection_uplift = prm.get_bool("Uplift and advect with fastscape");
          node_tolerance = prm.get_double("Node tolerance");
          noise_elevation = prm.get_double("Initial noise magnitude");
          sediment_rain_rates = Utilities::string_to_double
                                (Utilities::split_string_list(prm.get ("Sediment rain rates")));
          sediment_rain_times = Utilities::string_to_double
                                (Utilities::split_string_list(prm.get ("Sediment rain time intervals")));

          if (!this->convert_output_to_years())
            {
              maximum_fastscape_timestep /= year_in_seconds;
              for (unsigned int j=0; j<sediment_rain_rates.size(); ++j)
                sediment_rain_rates[j] *= year_in_seconds;
            }

          if (sediment_rain_rates.size() != sediment_rain_times.size()+1)
            AssertThrow(false, ExcMessage("Error: There must be one more sediment rain rate than time interval."));

          for (unsigned int i=1; i<sediment_rain_times.size(); ++i)
            AssertThrow(sediment_rain_times[i] > sediment_rain_times[i-1], ExcMessage("Sediment rain time intervals must be an increasing array."));

          prm.enter_subsection("Boundary conditions");
          {
            bottom = prm.get_integer("Front");
            right = prm.get_integer("Right");
            top = prm.get_integer("Back");
            left = prm.get_integer("Left");
            left_flux = prm.get_double("Left mass flux");
            right_flux = prm.get_double("Right mass flux");
            top_flux = prm.get_double("Back mass flux");
            bottom_flux = prm.get_double("Front mass flux");

            if (!this->convert_output_to_years())
              {
                left_flux *= year_in_seconds;
                right_flux *= year_in_seconds;
                top_flux *= year_in_seconds;
                bottom_flux *= year_in_seconds;
              }

            // Put the boundary condition values into a four digit value to send to FastScape.
            fastscape_boundary_conditions = bottom*1000+right*100+top*10+left;

            if ((left_flux != 0 && top_flux != 0) || (left_flux != 0 && bottom_flux != 0) ||
                (right_flux != 0 && bottom_flux != 0) || (right_flux != 0 && top_flux != 0))
              AssertThrow(false,ExcMessage("Currently the plugin does not support mass flux through adjacent boundaries."));

            topbottom_ghost_nodes_periodic = prm.get_bool("Back front ghost nodes periodic");
            leftright_ghost_nodes_periodic = prm.get_bool("Left right ghost nodes periodic");
          }
          prm.leave_subsection();

          prm.enter_subsection("Erosional parameters");
          {
            drainage_area_exponent_m = prm.get_double("Drainage area exponent");
            slope_exponent_n = prm.get_double("Slope exponent");
            sediment_river_incision_rate = prm.get_double("Sediment river incision rate");
            bedrock_river_incision_rate = prm.get_double("Bedrock river incision rate");
            sediment_transport_coefficient = prm.get_double("Sediment diffusivity");
            bedrock_transport_coefficient = prm.get_double("Bedrock diffusivity");
            bedrock_deposition_g = prm.get_double("Bedrock deposition coefficient");
            sediment_deposition_g = prm.get_double("Sediment deposition coefficient");
            slope_exponent_p = prm.get_double("Multi-direction slope exponent");
            flat_elevation = prm.get_integer("Orographic elevation control");
            wind_barrier_elevation = prm.get_integer("Orographic wind barrier height");
            flat_erosional_factor = prm.get_double("Elevation factor");
            wind_barrier_erosional_factor = prm.get_double("Wind barrier factor");
            stack_controls = prm.get_bool("Stack orographic controls");
            use_orographic_controls = prm.get_bool("Flag to use orographic controls");

            if (!this->convert_output_to_years())
              {
                bedrock_river_incision_rate *= year_in_seconds;
                bedrock_transport_coefficient *= year_in_seconds;
                sediment_river_incision_rate *= year_in_seconds;
                bedrock_transport_coefficient *= year_in_seconds;
              }

            // Wind direction
            if (prm.get ("Wind direction") == "west")
              wind_direction = 0;
            else if (prm.get ("Wind direction") == "east")
              wind_direction = 1;
            else if (prm.get ("Wind direction") == "north")
              wind_direction = 2;
            else if (prm.get ("Wind direction") == "south")
              wind_direction = 3;
            else
              AssertThrow(false, ExcMessage("Not a valid wind direction."));

            // set fixed ghost nodes to a base level for erosion that differs from sea level
            use_fixed_erosional_base = prm.get_bool("Use a fixed erosional base level");
            if (use_fixed_erosional_base)
              AssertThrow(use_fixed_erosional_base && use_ghost_nodes, ExcMessage(
                            "If you want to use an erosional base level differing from sea level, "
                            "you need to use ghost nodes."));
            h_erosional_base = prm.get_double("Erosional base level");
          }
          prm.leave_subsection();

          prm.enter_subsection("Marine parameters");
          {
            sea_level = prm.get_double("Sea level");
            sand_surface_porosity = prm.get_double("Sand porosity");
            silt_surface_porosity = prm.get_double("Silt porosity");
            sand_efold_depth = prm.get_double("Sand e-folding depth");
            silt_efold_depth = prm.get_double("Silt e-folding depth");
            sand_silt_ratio = prm.get_double("Sand-silt ratio");
            sand_silt_averaging_depth = prm.get_double("Depth averaging thickness");
            sand_transport_coefficient = prm.get_double("Sand transport coefficient");
            silt_transport_coefficient = prm.get_double("Silt transport coefficient");

            if (!this->convert_output_to_years())
              {
                sand_transport_coefficient *= year_in_seconds;
                silt_transport_coefficient *= year_in_seconds;
              }
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection ();

      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Visualization");
        {
          output_interval = prm.get_double ("Time between graphical output");
          if (this->convert_output_to_years())
            output_interval *= year_in_seconds;
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace MeshDeformation
  {
    ASPECT_REGISTER_MESH_DEFORMATION_MODEL(FastScape,
                                           "fastscape",
                                           "A plugin that uses the program FastScape to compute the deformation of the mesh surface. "
                                           "FastScape is a surface processes code that computes the erosion, transport and "
                                           "deposition of sediments both on land and in the marine domain. These surface processes "
                                           "include river incision (through the stream power law), hillslope diffusion and "
                                           "marine diffusion, as described in Braun and Willett 2013; Yuan et al. 2019; "
                                           "Yuan et al. 2019b. "
                                           "\n"
                                           "Upon initialization, FastScape requires the initial topography of the surface boundary "
                                           "of ASPECT's model domain and several user-specified erosional and depositional parameters. "
                                           "In each ASPECT timestep, FastScape is then fed ASPECT's material velocity at the surface "
                                           "boundary. The z-component of this velocity is used to uplift the FastScape surface, "
                                           "while the horizontal components are used to advect the topography in the x-y plane. "
                                           "\n"
                                           "After solving its governing equations (this can be done in several timesteps "
                                           "that are smaller than the ASPECT timestep), FastScape returns a new topography of the surface. "
                                           "The difference in topography before and after the call to FastScape divided by the ASPECT "
                                           "timestep provides the mesh velocity at the domain's surface that is used to displace the surface "
                                           "and internal mesh nodes. "
                                           "\n"
                                           "FastScape can be used in both 2D and 3D ASPECT simulations. In 2D, one can think of the coupled "
                                           "model as a T-model.The ASPECT domain spans the x - z plane, while FastScape acts on the horizontal "
                                           "x-y plane. This means that to communicate ASPECT's material velocities to FastScape, "
                                           "FastScape mesh nodes with the same x-coordinate (so lying along the y-direction) get the same velocities. "
                                           "In turn, the FastScape topography is collapsed back onto the line of the ASPECT surface boundary "
                                           "by averaging the topography over the y-direction. In 3D no such actions are necessary. "
                                           "\n"
                                           "The FastScape manual (https://fastscape.org/fastscapelib-fortran/) provides more information "
                                           "on the input parameters. ")

  }
}
#endif
