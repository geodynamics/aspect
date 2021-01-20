/*
  Copyright (C) 2019 by the authors of the ASPECT code.
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

#include <aspect/mesh_deformation/fastscape.h>
#include <aspect/geometry_model/box.h>
#include <deal.II/numerics/vector_tools.h>
#include <ctime>

namespace aspect
{
  namespace MeshDeformation
  {

    template <int dim>
    void
    FastScape<dim>::initialize ()
    {
      AssertThrow(Plugins::plugin_type_matches<const GeometryModel::Box<dim>>(this->get_geometry_model()),
                  ExcMessage("Fastscape can only be run with a box model"));

      const GeometryModel::Box<dim> *geometry
        = dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model());

      // Find the id associated with the top boundary and boundaries that call mesh deformation.
      const types::boundary_id relevant_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id ("top");
      const std::set<types::boundary_id> boundary_ids
        = this->get_mesh_deformation_handler().get_active_mesh_deformation_boundary_indicators();

      // Get the deformation type names called for each boundary.
      /*
       * TODO: Why does this only work if I declare it in the header file?
       * If I declare it here, then I get a no match operand[] error when trying to get the
       * names variable.
       */
      mesh_deformation_boundary_indicators_map
        = this->get_mesh_deformation_handler().get_active_mesh_deformation_names();

      // Loop over each mesh deformation boundary, and make sure FastScape is only called on the surface.
      for (std::set<types::boundary_id>::const_iterator p = boundary_ids.begin();
           p != boundary_ids.end(); ++p)
        {
          const std::vector<std::string> names = mesh_deformation_boundary_indicators_map[*p];
          for (unsigned int i = 0; i < names.size(); ++i )
            {
              if (names[i] == "fastscape")
                AssertThrow((*p == relevant_boundary),
                            ExcMessage("Fastscape can only be called on the surface boundary."));
            }
        }

      // Initialize parameters for restarting fastscape
      restart = this->get_parameters().resume_computation;
      restart_step = 0;

      // second is for maximum coordinates, first for minimum.
      for (unsigned int i=0; i<dim; ++i)
        {
          grid_extent[i].first = geometry->get_origin()[i];
          grid_extent[i].second = geometry->get_extents()[i];
        }

      // TODO: There has to be a better type to use to get this.
      // const unsigned int repetitions[dim] = {geometry->get_repetitions()};
      const unsigned int x_repetitions = geometry->get_repetitions(0);
      const unsigned int y_repetitions = geometry->get_repetitions(1);

      // Set nx and dx, as these will be the same regardless of dimension.
      nx = 1+(2*use_ghost)+std::pow(2,surface_resolution+additional_refinement)*x_repetitions;
      dx = (grid_extent[0].second - grid_extent[0].first)/(nx-1-(2*use_ghost));
      x_extent = (grid_extent[0].second - grid_extent[0].first)+2*dx*use_ghost;

      // Sub intervals are 3 less than points, if including the ghost nodes. Otherwise 1 less.
      table_intervals[0] = nx-1-(2*use_ghost);
      // TODO: it'd be best to not have to use dim-1 intervals at all.
      table_intervals[dim-1] = 1;

      if (dim == 2)
        {
          dy = dx;
          y_extent = round(y_extent_2d/dy)*dy+2*dy*use_ghost;
          ny = 1+y_extent/dy;
        }
      else
        {
          ny = 1+(2*use_ghost)+std::pow(2,surface_resolution+additional_refinement)*y_repetitions;
          dy = (grid_extent[1].second - grid_extent[1].first)/(ny-1-(2*use_ghost));
          table_intervals[1] = ny-1-(2*use_ghost);
          y_extent = (grid_extent[1].second - grid_extent[1].first)+2*dy*use_ghost;
        }

      // Determine array size to send to fastscape
      array_size = nx*ny;

      // Create a folder for the FastScape visualization files.
      Utilities::create_directory (this->get_output_directory() + "VTK/",
                                   this->get_mpi_communicator(),
                                   false);
    }


    template <int dim>
    void
    FastScape<dim>::compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                             ConstraintMatrix &mesh_velocity_constraints,
                                                             const std::set<types::boundary_id> &boundary_ids) const
    {
      TimerOutput::Scope timer_section(this->get_computing_timer(), "Fastscape plugin");
      const types::boundary_id relevant_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id ("top");
      const int current_timestep = this->get_timestep_number ();

      double a_dt = this->get_timestep();
      if (this->convert_output_to_years())
        {
          a_dt = this->get_timestep()/year_in_seconds;
        }

      // We only want to run fastscape if there was a change in time.
      if (a_dt > 0)
        {
          /*
           * Initialize a vector of temporary variables to hold: z component, index, Vx, Vy, and Vz.
           */
          std::vector<std::vector<double>> temporary_variables(dim+2, std::vector<double>());
          //Vector to hold the velocities that represent the change to the surface.
          std::vector<double> V(array_size);
          //precision = 0.001;

          // Get a quadrature rule that exists only on the corners, and increase the refinement if specified.
          const QIterated<dim-1> face_corners (QTrapez<1>(),
                                               pow(2,additional_refinement+resolution_difference));

          FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                            this->get_fe(),
                                            face_corners,
                                            update_values |
                                            update_quadrature_points);

          typename DoFHandler<dim>::active_cell_iterator
          cell = this->get_dof_handler().begin_active(),
          endc = this->get_dof_handler().end();

          for (; cell != endc; ++cell)
            if (cell->is_locally_owned() && cell->at_boundary())
              for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
                if (cell->face(face_no)->at_boundary())
                  {
                    if ( cell->face(face_no)->boundary_id() != relevant_boundary)
                      continue;

                    std::vector<Tensor<1,dim> > vel(face_corners.size());
                    fe_face_values.reinit(cell, face_no);
                    fe_face_values[this->introspection().extractors.velocities].get_function_values(this->get_solution(), vel);

                    for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
                      {
                        const Point<dim> vertex = fe_face_values.quadrature_point(corner);

                        // Find what x point we're at. Add 1 or 2 depending on if ghost nodes are used.
                        const double indx = 1+use_ghost+vertex(0)/dx;

                        // If our x or y index isn't close to a whole number, then it's likely an artifact
                        // from using an over-resolved quadrature rule, in that case ignore it.
                        if (abs(indx - round(indx)) >= precision)
                          continue;

                        /*
                         * If we're in 2D, we want to take the values and apply them to every row of X points.
                         */
                        if (dim == 2)
                          {
                            for (int ys=0; ys<ny; ys++)
                              {
                                /*
                                 * Fastscape indexes from 1 to n, starting at X and Y = 0, and increases
                                 * across the X row. At the end of the row, it jumps back to X = 0
                                 * and up to the next X row in increasing Y direction. We track
                                 * this to correctly place the variables later on.
                                 * Nx*ys effectively tells us what row we are in
                                 * and then indx tells us what position in that row.
                                 */
                                const double index = round(indx)+nx*ys;

                                temporary_variables[0].push_back(vertex(dim-1) - grid_extent[dim-1].second);
                                temporary_variables[1].push_back(index-1);

                                for (unsigned int i=0; i<dim; ++i)
                                  {
                                    if (this->convert_output_to_years())
                                      temporary_variables[i+2].push_back(vel[corner][i]*year_in_seconds);
                                    else
                                      temporary_variables[i+2].push_back(vel[corner][i]);
                                  }
                              }
                          }
                        else
                          {
                            //Because indy only gives us the row we're in, we don't need to add 2 for the ghost node.
                            const double indy = 1+use_ghost+vertex(1)/dy;
                            if (abs(indy - round(indy)) >= precision)
                              continue;

                            const double index = round((indy-1))*nx+round(indx);

                            temporary_variables[0].push_back(vertex(dim-1) - grid_extent[dim-1].second);   //z component
                            temporary_variables[1].push_back(index-1);

                            for (unsigned int i=0; i<dim; ++i)
                              {
                                if (this->convert_output_to_years())
                                  temporary_variables[i+2].push_back(vel[corner][i]*year_in_seconds);
                                else
                                  temporary_variables[i+2].push_back(vel[corner][i]);
                              }
                          }
                      }
                  }

          // Run fastscape on single processor.
          if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            {
              /*
               * Initialize the variables that will be sent to fastscape.
               * These have to be doubles of array_size, which C++ doesn't like,
               * so they're initialized this way.
               */
              std::unique_ptr<double[]> h (new double[array_size]());
              std::unique_ptr<double[]> vx (new double[array_size]());
              std::unique_ptr<double[]> vy (new double[array_size]());
              std::unique_ptr<double[]> vz (new double[array_size]());
              std::unique_ptr<double[]> kf (new double[array_size]());
              std::unique_ptr<double[]> kd (new double[array_size]());
              std::unique_ptr<double[]> slopep (new double[array_size]());
              std::unique_ptr<double[]> b (new double[array_size]());
              std::vector<double> h_old(array_size);


              // Create variables for output directory and restart file
              std::string dirname = this->get_output_directory();
              const char *c=dirname.c_str();
              int length = dirname.length();
              const std::string restart_filename = dirname + "fastscape_h_restart.txt";
              const std::string restart_step_filename = dirname + "fastscape_steps_restart.txt";
              const std::string restart_filename_basement = dirname + "fastscape_b_restart.txt";

              // Initialize kf and kd.
              for (int i=0; i<array_size; i++)
                {
                  kf[i] = kff;
                  kd[i] = kdd;
                }

              // Get info from first processor.
              for (unsigned int i=0; i<temporary_variables[1].size(); i++)
                {
                  h[temporary_variables[1][i]]= temporary_variables[0][i];
                  vx[temporary_variables[1][i]]= temporary_variables[2][i];
                  vz[temporary_variables[1][i]]= temporary_variables[dim+1][i];

                  if (dim == 2 )
                    vy[temporary_variables[1][i]]=0;
                  else
                    vy[temporary_variables[1][i]]=temporary_variables[3][i];
                }

              for (unsigned int p=1; p<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++p)
                {
                  // First, find out the size of the array a processor wants to send.
                  MPI_Status status;
                  MPI_Probe(p, 42, this->get_mpi_communicator(), &status);
                  int incoming_size = 0;
                  MPI_Get_count(&status, MPI_DOUBLE, &incoming_size);

                  // Resize the array so it fits whatever the processor sends.
                  for (unsigned int i=0; i<temporary_variables.size(); ++i)
                    {
                      temporary_variables[i].resize(incoming_size);
                    }

                  for (unsigned int i=0; i<temporary_variables.size(); i++)
                    MPI_Recv(&temporary_variables[i][0], incoming_size, MPI_DOUBLE, p, 42, this->get_mpi_communicator(), &status);

                  // Now, place the numbers into the correct place based off the index.
                  for (unsigned int i=0; i<temporary_variables[1].size(); i++)
                    {
                      h[temporary_variables[1][i]] = temporary_variables[0][i];
                      vx[temporary_variables[1][i]] = temporary_variables[2][i];
                      vz[temporary_variables[1][i]] = temporary_variables[dim+1][i];

                      //In a 2D case, we don't know any y velocities.
                      if (dim == 2 )
                        vy[temporary_variables[1][i]] = 0;
                      else
                        vy[temporary_variables[1][i]] = temporary_variables[3][i];
                    }
                }

              if (current_timestep == 1 || restart)
                {
                  this->get_pcout() << "   Initializing FastScape... " << (1+surface_resolution+additional_refinement) <<
                                    " levels, cell size: " << dx << " m." << std::endl;

                  // If we are restarting from a checkpoint, load h values for fastscape so we don't lose resolution.
                  if (restart)
                    {
                      this->get_pcout() << "      Loading FastScape restart file... " << std::endl;

                      // Load in h values.
                      std::ifstream in;
                      in.open(restart_filename.c_str());
                      if (in)
                        {
                          int line = 0;

                          while (line < array_size)
                            {
                              in >> h[line];
                              line++;
                            }

                          in.close();
                        }
                      else if (!in)
                        AssertThrow(false,ExcMessage("Cannot open file to restart FastScape."));

                      // Load in b values.
                      std::ifstream in_b;
                      in_b.open(restart_filename_basement.c_str());
                      if (in_b)
                        {
                          int line = 0;

                          while (line < array_size)
                            {
                              in_b >> b[line];
                              line++;
                            }

                          in_b.close();
                        }
                      else if (!in_b)
                        AssertThrow(false,ExcMessage("Cannot open basement file to restart FastScape."));

                      /*
                       * Now load the fastscape istep at time of restart.
                       * Reinitializing fastscape always resets this to 0, so here
                       * we keep it in a separate variable to keep track for visualization files.
                       */
                      std::ifstream in_step;
                      in_step.open(restart_step_filename.c_str());
                      if (in_step)
                        {
                          in_step >> restart_step;
                          in_step.close();
                        }
                      else if (!in_step)
                        AssertThrow(false,ExcMessage("Cannot open file to restart fastscape."));

                      restart = false;
                    }

                  // Initialize fastscape with grid and extent.
                  fastscape_init_();
                  fastscape_set_nx_ny_(&nx,&ny);
                  fastscape_setup_();
                  fastscape_set_xl_yl_(&x_extent,&y_extent);

                  // Set boundary conditions
                  fastscape_set_bc_(&bc);

                  // Initialize topography
                  fastscape_init_h_(h.get());

                  // Set erosional parameters. May have to move this if sed values are updated over time.
                  fastscape_set_erosional_parameters_(kf.get(), &kfsed, &m, &n, kd.get(), &kdsed, &g, &gsed, &p);

                  if (use_marine)
                    fastscape_set_marine_parameters_(&sl, &p1, &p2, &z1, &z2, &r, &l, &kds1, &kds2);

                  // Only set the basement if it's a restart
                  if (current_timestep != 1)
                    fastscape_set_basement_(b.get());
                }
              else
                {
                  // If it isn't the first timestep we ignore initialization and instead copy all height values from FastScape.
                 if(use_v)
                   fastscape_copy_h_(h.get());
                }

              /*
               * Generally, any additional modifications should occur below this point.
               */

              /*
               * The ghost nodes are added as a single layer of nodes surrounding the entire model,
               * so that 3x3 amount of nodes representing ASPECT will become 5x5 in FastScape, with the inner
               * 3x3 being returned to ASPECT.
               */
              // I redid the indexing here, at some point I should double check that these still work without issue.
              if (use_ghost)
                {
                  /*
                  * Copy the slopes at each point, this will be used to set an H
                  * at the ghost nodes if a boundary mass flux is given.
                  */
                  fastscape_copy_slope_(slopep.get());

                  /*
                   * Here we set the ghost nodes at the left and right boundaries. In most cases,
                   * this involves setting the node to the same values of v and h as the inward node.
                   * With the inward node being above or below for the bottom and top rows of ghost nodes,
                   * or to the left and right for the right and left columns of ghost nodes.
                   */
                  for (int j=0; j<ny; j++)
                    {
                      /*
                      * Nx*j will give us the row we're in, and one is subtracted as FastScape starts from 1 not zero.
                      * If we're on the left, the multiple of the row will always represent the first node.
                      * Subtracting one to the row above this gives us the last node of the previous row.
                      */
                      const int index_left = nx*j;
                      const int index_right = nx*(j+1)-1;
                      double slope = 0;

                      /*
                      * Here we set the ghost nodes to the ones next to them, where for the left we
                      * add one to go to the node to the right, and for the right side
                      * we subtract one to go to the inner node to the left.
                      */
                      vz[index_right] = vz[index_right-1];
                      vz[index_left] =  vz[index_left+1];

                      vy[index_right] = vy[index_right-1];
                      vy[index_left] = vy[index_left+1];

                      vx[index_right] = vx[index_right-1];
                      vx[index_left] = vx[index_left+1];

                      if (current_timestep == 1 || left_flux == 0)
                        {
                          /*
                           * If its the first timestep add in initial slope. If we have no flux,
                           * set the ghost node to the node next to it.
                           * FastScape calculates the slope by looking at all nodes surrounding the point
                           * so we need to consider the slope over 2 dx.
                           */
                          slope = left_flux/kdd;
                          h[index_left] = h[index_left+1] + slope*2*dx;
                        }
                      else
                        {
                          /*
                           * If we have flux through a boundary, we need to update the height to keep the correct slope.
                           * Because the corner nodes always show a slope of zero, this will update them according to
                           * the closest non-ghost node. E.g. if we're at a corner node, look instead up a row and inward.
                           */
                          if (j == 0)
                            slope = left_flux/kdd - std::tan(slopep[index_left+nx+1]*numbers::PI/180.);
                          else if (j==(ny-1))
                            slope = left_flux/kdd - std::tan(slopep[index_left-nx+1]*numbers::PI/180.);
                          else
                            slope = left_flux/kdd - std::tan(slopep[index_left+1]*numbers::PI/180.);

                          h[index_left] = h[index_left] + slope*2*dx;
                        }

                      if (current_timestep == 1 || right_flux == 0)
                        {
                          slope = right_flux/kdd;
                          h[index_right] = h[index_right-1] + slope*2*dx;
                        }
                      else
                        {
                          if (j == 0)
                            slope = right_flux/kdd - std::tan(slopep[index_right+nx-1]*numbers::PI/180.);
                          else if (j==(ny-1))
                            slope = right_flux/kdd - std::tan(slopep[index_right-nx-1]*numbers::PI/180.);
                          else
                            slope = right_flux/kdd - std::tan(slopep[index_right-1]*numbers::PI/180.);

                          h[index_right] = h[index_right] + slope*2*dx;
                        }

                      /*
                      * If the boundaries are periodic, then we look at the velocities on both sides of the
                      * model, and set the ghost node according to the direction of flow. As FastScape will
                      * receive all velocities it will have a direction, and we only need to look at the (non-ghost)
                      * nodes directly to the left and right.
                      */
                      if (left == 0 && right == 0)
                        {
                          // First we assume that flow is going to the left.
                          int side = index_left;

                          // Indexing depending on which side the ghost node is being set to.
                          int jj = 1;

                          /*
                          * If nodes on both sides are going the same direction, then set the respective
                          * ghost nodes to equal these sides. By doing this, the ghost nodes at the opposite
                          * side of flow will work as a mirror mimicing what is happening on the other side.
                          */
                          if (vx[index_right-1] > 0 && vx[index_left+1] >= 0)
                            {
                              side = index_right;
                              jj = -1;
                            }
                          else if (vx[index_right-1] <= 0 && vx[index_left+1] < 0)
                            {
                              side = index_left;
                              jj = 1;
                            }
                          else
                            continue;

                          vz[index_right] = vz[side+jj];
                          vz[index_left] = vz[side+jj];

                          vy[index_right] = vy[side+jj];
                          vy[index_left] = vy[side+jj];

                          vx[index_right] = vx[side+jj];
                          vx[index_left] = vx[side+jj];

                          h[index_right] = h[side+jj];
                          h[index_left] = h[side+jj];
                        }
                    }

                  // Now do the same for the top and bottom ghost nodes.
                  for (int j=0; j<nx; j++)
                    {
                      // The bottom row indexes are 0 to nx-1.
                      const int index_bot = j;

                      // Nx multiplied by (total rows - 1) gives us the start of
                      // the top row, and j gives the position in the row.
                      const int index_top = nx*(ny-1)+j;
                      double slope = 0;

                      vz[index_bot] = vz[index_bot+nx];
                      vz[index_top] = vz[index_top-nx];

                      vy[index_bot] = vy[index_bot+nx];
                      vy[index_top] =  vy[index_top-nx];

                      vx[index_bot] = vx[index_bot+nx];
                      vx[index_top] =  vx[index_top-nx];

                      if (current_timestep == 1 || top_flux == 0)
                        {
                          slope = top_flux/kdd;
                          h[index_top] = h[index_top-nx] + slope*2*dx;
                        }
                      else
                        {
                          if (j == 0)
                            slope = top_flux/kdd - std::tan(slopep[index_top-nx+1]*numbers::PI/180.);
                          else if (j==(nx-1))
                            slope = top_flux/kdd - std::tan(slopep[index_top-nx-1]*numbers::PI/180.);
                          else
                            slope = top_flux/kdd - std::tan(slopep[index_top-nx]*numbers::PI/180.);

                          h[index_top] = h[index_top] + slope*2*dx;
                        }

                      if (current_timestep == 1 || bottom_flux == 0)
                        {
                          slope = bottom_flux/kdd;
                          h[index_bot] = h[index_bot+nx] + slope*2*dx;
                        }
                      else
                        {
                          if (j == 0)
                            slope = bottom_flux/kdd - std::tan(slopep[index_bot+nx+1]*numbers::PI/180.);
                          else if (j==(nx-1))
                            slope = bottom_flux/kdd - std::tan(slopep[index_bot+nx-1]*numbers::PI/180.);
                          else
                            slope = bottom_flux/kdd - std::tan(slopep[index_bot+nx]*numbers::PI/180.);

                          h[index_bot] = h[index_bot] + slope*2*dx;
                        }

                      if (bottom == 0 && top == 0)
                        {
                          int side = index_bot;
                          int jj = nx;

                          if (vy[index_bot+nx-1] > 0 && vy[index_top-nx-1] >= 0)
                            {
                              side = index_top;
                              jj = -nx;
                            }
                          else if (vy[index_bot+nx-1] <= 0 && vy[index_top-nx-1] < 0)
                            {
                              side = index_bot;
                              jj = nx;
                            }
                          else
                            continue;

                          vz[index_bot] = vz[side+jj];
                          vz[index_top] = vz[side+jj];

                          vy[index_bot] = vy[side+jj];
                          vy[index_top] =  vy[side+jj];

                          vx[index_bot] = vx[side+jj];
                          vx[index_top] =  vx[side+jj];

                          h[index_bot] = h[side+jj];
                          h[index_top] = h[side+jj];
                        }
                    }
                }


                // Find the appropriate sediment rain based off the time interval.
                double time = this->get_time()/year_in_seconds;
                double sediment_rain = sr_values[0];
                for (unsigned int j=0; j<sr_times.size(); j++)
                {
                  if(time > sr_times[j])
                       sediment_rain = sr_values[j+1];
                }

              /*
               * Keep initial h values so we can calculate velocity later.
               * In the first timestep, h will be given from other processors.
               * In later timesteps, we copy h directly from fastscape.
               */
              std::srand(fs_seed);
              for (int i=0; i<array_size; i++)
                {
                  h_old[i] = h[i];

                  // Initialize random noise after h_old is set, so aspect sees this initial topography change.
                  if (current_timestep == 1)
                    {
                      // + or - 5 meters of topography.
                      const double h_seed = (std::rand()%100)/10 - 5;
                      h[i] = h[i] + h_seed;
                    }
              
                  // Here we add the sediment rain (m/yr) as a flat increase in height.
                  // This is done because adding it as an uplift rate would affect the basement.
                  if (sediment_rain > 0)
                    {
                      // Only apply sediment rain to areas below sea level.
                      if (h[i] < sl)
                        {
                          // If the rain would put us above sea level, set height to sea level.
                          if (h[i] + sediment_rain*a_dt > sl)
                            h[i] = sl;
                          else
                            h[i] = h[i] + sediment_rain*a_dt;
                        }
                    }
                }

              // Get current fastscape timestep.
              int istep = 0;
              fastscape_get_step_(&istep);

              // Write a file to store h & step in case of restart.
              // TODO: there's probably a faster way to write these.
              if ((this->get_parameters().checkpoint_time_secs == 0) &&
                  (this->get_parameters().checkpoint_steps > 0) &&
                  (current_timestep % this->get_parameters().checkpoint_steps == 0))
                {

                  std::ofstream out_h (restart_filename.c_str());
                  std::ofstream out_step (restart_step_filename.c_str());
                  std::ofstream out_b (restart_filename_basement.c_str());
                  std::stringstream bufferb;
                  std::stringstream bufferh;

                  fastscape_copy_basement_(b.get());

                  out_step << (istep+restart_step) << "\n";

                  for (int i=0; i<array_size; i++)
                    {
                      bufferh << h[i] << "\n";
                      bufferb << b[i] << "\n";
                    }

                  out_h << bufferh.str();
                  out_b << bufferb.str();
                }

              // Find a fastscape timestep that is below our maximum timestep.
              int steps = nstep;
              double f_dt = a_dt/steps;
              while (f_dt>max_timestep)
                {
                  steps=steps*2;
                  f_dt = a_dt/steps;
                }

              // Set time step
              fastscape_set_dt_(&f_dt);

              // Set velocity components and h.
              if(use_v)
              {
                fastscape_set_u_(vz.get());
                fastscape_set_v_(vx.get(), vy.get());
              }

              fastscape_set_h_(h.get());

              // The visualization number depends on the FastScape step, and since this goes back
              // to zero after a restart, we need to keep it up to date.
              // TODO: it may be better to instead have visualization track the time.
              int visualization_step = istep+restart_step;
              steps = istep+steps;

              this->get_pcout() << "   Calling FastScape... " << (steps-istep) << " timesteps of " << f_dt << " years." << std::endl;
              {
                auto t_start = std::chrono::high_resolution_clock::now();

                /*
                 * If we use stratigraphy it'll handle visualization and not the normal function.
                 * TODO: The frequency in this needs to be the same as the total timesteps fastscape will
                 * run for, need to figure out how to work this in better.
                 */
                if (use_strat && current_timestep == 1)
                  fastscape_strati_(&nstepp, &nreflectorp, &steps, &vexp);
                else if (!use_strat)
                  fastscape_named_vtk_(h.get(), &vexp, &visualization_step, c, &length);

                do
                  {
                    // Execute step, this increases timestep counter
                    fastscape_execute_step_();

                    // Get value of time step counter
                    fastscape_get_step_(&istep);

                    // Outputs new h values
                    fastscape_copy_h_(h.get());
                  }
                while (istep<steps);

                // Output out how long FastScape took to run.
                auto t_end = std::chrono::high_resolution_clock::now();
                double r_time = std::chrono::duration<double>(t_end-t_start).count();
                this->get_pcout() << "      FastScape runtime... " << round(r_time*1000)/1000 << "s" << std::endl;
              }

              // If we've reached the end time, destroy fastscape.
              if (this->get_time()+a_dt >= end_time)
                {
                  this->get_pcout() << "      Destroying FastScape..." << std::endl;

                  // We generally create the visualization file before running FastScape,
                  // because we won't run it again we want to output the final result.
                  visualization_step = visualization_step+1;
                  if (!use_strat)
                    fastscape_named_vtk_(h.get(), &vexp, &visualization_step, c, &length);

                  fastscape_destroy_();
                }

              // Find out our velocities from the change in height.
              // Where V is a vector of array size that exists on all processors.
              for (int i=0; i<array_size; i++)
                {
                  V[i] = (h[i] - h_old[i])/a_dt;
                }

              MPI_Bcast(&V[0], array_size, MPI_DOUBLE, 0, this->get_mpi_communicator());
            }
          else
            {
              for (unsigned int i=0; i<temporary_variables.size(); i++)
                MPI_Ssend(&temporary_variables[i][0], temporary_variables[1].size(), MPI_DOUBLE, 0, 42, this->get_mpi_communicator());

              MPI_Bcast(&V[0], array_size, MPI_DOUBLE, 0, this->get_mpi_communicator());
            }

          // Get the sizes needed for the data table.
          TableIndices<dim> size_idx;
          for (unsigned int d=0; d<dim; ++d)
            {
              size_idx[d] = table_intervals[d]+1;
            }

          // Initialize a table to hold all velocity values that will be interpolated back to ASPECT.
          Table<dim,double> data_table;
          data_table.TableBase<dim,double>::reinit(size_idx);
          TableIndices<dim> idx;

          // Loop through the data table and fill out the surface (i or k = 1) to the velocities from FastScape.
          if (dim == 2)
            {
              std::vector<double> V2(nx);

              for (int i=1; i<(nx-1); i++)
                {
                  // If using the center slice, find velocities from the row closest to the center.
                  if (slice)
                    {
                      const int index = i+nx*(round((ny-1)/2));
                      V2[i-1] = V[index];
                    }
                  // Here we use average velocities across the y nodes, excluding the ghost nodes (top and bottom row).
                  // Note: If ghost nodes are turned off, boundary effects may influence this.
                  else
                    {
                      for (int ys=1; ys<(ny-1); ys++)
                        {
                          const int index = i+nx*ys;
                          V2[i-1] += V[index];
                        }
                      V2[i-1] = V2[i-1]/(ny-2);
                    }
                }

              for (unsigned int i=0; i<data_table.size()[1]; ++i)
                {
                  idx[1] = i;

                  for (unsigned int j=0; j<(data_table.size()[0]); ++j)
                    {
                      idx[0] = j;

                      // Convert from m/yr to m/s if needed. i == 1 is the surface and i = 0 the bottom.
                      if (i == 1)
                        {
                          if (this->convert_output_to_years())
                            data_table(idx) = V2[j]/year_in_seconds;
                          else
                            data_table(idx) = V2[j];
                        }
                      else
                        data_table(idx)= 0;
                    }
                }
            }
          else
            {
              // Indexes through z, y, and then x.
              for (unsigned int k=0; k<data_table.size()[2]; ++k)
                {
                  idx[2] = k;

                  for (unsigned int i=0; i<data_table.size()[1]; ++i)
                    {
                      idx[1] = i;

                      for (unsigned int j=0; j<data_table.size()[0]; ++j)
                        {
                          idx[0] = j;

                          // Tables are initialized with k = 0 as the bottom and 1 as top, we are only interested in the surface so only fill k==1.
                          if (k==1)
                            {
                              // Fill table, where nx+1 allows skipping of the first ghost node.
                              if (this->convert_output_to_years())
                                data_table(idx) = V[(nx+1)*use_ghost+nx*i+j]/year_in_seconds;
                              else
                                data_table(idx) = V[(nx+1)*use_ghost+nx*i+j];
                            }
                          else
                            data_table(idx)= 0;
                        }
                    }
                }
            }

          Functions::InterpolatedUniformGridData<dim> *velocities;
          velocities = new Functions::InterpolatedUniformGridData<dim> (grid_extent,
                                                                        table_intervals,
                                                                        data_table);

          auto lambda = [&](const Point<dim> &p) -> double
          {
            return velocities->value(p);
          };

          VectorFunctionFromScalarFunctionObject<dim> vector_function_object(
            lambda,
            dim-1,
            dim);

          VectorTools::interpolate_boundary_values (mesh_deformation_dof_handler,
                                                    *boundary_ids.begin(),
                                                    vector_function_object,
                                                    mesh_velocity_constraints);
        }
    }


    // TODO: Give better explanations of variables and cite the fastscape documentation.
    template <int dim>
    void FastScape<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection ("Fastscape");
        {
          prm.declare_entry("Number of steps", "10",
                            Patterns::Integer(),
                            "Number of steps per ASPECT timestep");
          prm.declare_entry("Maximum timestep", "10e3",
                            Patterns::Double(0),
                            "Maximum timestep for FastScape.");
          prm.declare_entry("Vertical exaggeration", "-1",
                            Patterns::Double(),
                            "Vertical exaggeration for FastScape's VTK file. -1 outputs topography, basement, and sealevel.");
          prm.declare_entry("Additional fastscape refinement", "0",
                            Patterns::Integer(),
                            "How many levels above ASPECT FastScape should be refined.");
          prm.declare_entry ("Use center slice for 2d", "false",
                             Patterns::Bool (),
                             "If this is set to true, then a 2D model will only consider the "
                             "center slice FastScape gives. If set to false, then aspect will"
                             "average the inner third of what FastScape calculates.");
          prm.declare_entry("Fastscape seed", "1000",
                            Patterns::Integer(),
                            "Seed used for adding an initial 0-20 m noise to FastScape topography.");
          prm.declare_entry("Surface resolution", "1",
                            Patterns::Integer(),
                            "This should be set to the highest ASPECT resolution level you expect at the surface.");
          prm.declare_entry("Resolution difference", "0",
                            Patterns::Integer(),
                            "The difference between the lowest and highest resolution level at surface. So if three resolution "
                            "levels are expected, this would be set to 2.");
          prm.declare_entry ("Use marine component", "false",
                             Patterns::Bool (),
                             "Flag to use the marine component of FastScape.");
          prm.declare_entry ("Use stratigraphy", "false",
                             Patterns::Bool (),
                             "Flag to use stratigraphy");
          prm.declare_entry("Total steps", "100000",
                            Patterns::Integer(),
                            "Total number of steps you expect in the FastScape model, only used if stratigraphy is turned on.");
          prm.declare_entry("Number of horizons", "1",
                            Patterns::Integer(),
                            "Number of horizons to track and visualize in FastScape.");
          prm.declare_entry("Y extent in 2d", "100000",
                            Patterns::Double(),
                            "Y extent when used with 2D");
          prm.declare_entry ("Use ghost nodes", "true",
                             Patterns::Bool (),
                             "Flag to use ghost nodes");
          prm.declare_entry ("Use velocities", "true",
                             Patterns::Bool (),
                             "Flag to use FastScape advection and uplift.");
          prm.declare_entry("Precision", "0.001",
                            Patterns::Double(),
                            "How close an ASPECT node needs to be to a FastScape node.");
          prm.declare_entry ("Sediment rain", "0",
                             Patterns::List (Patterns::Double(0)),
                             "Sediment rain values given as a list equal to the number of intervals.");
          prm.declare_entry ("Sediment rain intervals", "0",
                             Patterns::List (Patterns::Double(0)),
                             "A list of sediment rain times.");


          prm.enter_subsection ("Boundary conditions");
          {
            prm.declare_entry ("Bottom", "1",
                               Patterns::Integer (0, 1),
                               "Bottom boundary condition, where 1 is fixed and 0 is reflective.");
            prm.declare_entry ("Right", "1",
                               Patterns::Integer (0, 1),
                               "Right boundary condition, where 1 is fixed and 0 is reflective.");
            prm.declare_entry ("Top", "1",
                               Patterns::Integer (0, 1),
                               "Top boundary condition, where 1 is fixed and 0 is reflective.");
            prm.declare_entry ("Left", "1",
                               Patterns::Integer (0, 1),
                               "Left boundary condition, where 1 is fixed and 0 is reflective.");
            prm.declare_entry("Left mass flux", "0",
                              Patterns::Double(),
                              "Flux per unit length through left boundary (m^2/yr)");
            prm.declare_entry("Right mass flux", "0",
                              Patterns::Double(),
                              "Flux per unit length through right boundary (m^2/yr)");
            prm.declare_entry("Top mass flux", "0",
                              Patterns::Double(),
                              "Flux per unit length through top boundary (m^2/yr)");
            prm.declare_entry("Bottom mass flux", "0",
                              Patterns::Double(),
                              "Flux per unit length through bottom boundary (m^2/yr)");
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
            prm.declare_entry("Bedrock deposition coefficient", "-1",
                              Patterns::Double(),
                              "Deposition coefficient for bedrock");
            prm.declare_entry("Sediment deposition coefficient", "-1",
                              Patterns::Double(),
                              "Deposition coefficient for sediment");
            prm.declare_entry("Bedrock river incision rate", "-1",
                              Patterns::Double(),
                              "River incision rate for bedrock");
            prm.declare_entry("Sediment river incision rate", "-1",
                              Patterns::Double(),
                              "River incision rate for sediment ");
            prm.declare_entry("Bedrock diffusivity", "1",
                              Patterns::Double(),
                              "Diffusivity of bedrock.");
            prm.declare_entry("Sediment diffusivity", "-1",
                              Patterns::Double(),
                              "Diffusivity of sediment.");
          }
          prm.leave_subsection();

          prm.enter_subsection ("Marine parameters");
          {
            prm.declare_entry("Sea level", "0",
                              Patterns::Double(),
                              "Sea level in meters. ");
            prm.declare_entry("Sand porosity", "0.0",
                              Patterns::Double(),
                              "Porosity of sand. ");
            prm.declare_entry("Shale porosity", "0.0",
                              Patterns::Double(),
                              "Porosity of shale. ");
            prm.declare_entry("Sand e-folding depth", "1e3",
                              Patterns::Double(),
                              "E-folding depth for the exponential of the sand porosity law.");
            prm.declare_entry("Shale e-folding depth", "1e3",
                              Patterns::Double(),
                              "E-folding depth for the exponential of the shale porosity law.");
            prm.declare_entry("Sand-shale ratio", "0.5",
                              Patterns::Double(),
                              "Ratio of sand to shale for material leaving continent.");
            prm.declare_entry("Depth averaging thickness", "1e2",
                              Patterns::Double(),
                              "Depth averaging for the sand-shale equation in meters.");
            prm.declare_entry("Sand transport coefficient", "5e2",
                              Patterns::Double(),
                              "Transport coefficient for sand in m^2/yr ");
            prm.declare_entry("Shale transport coefficient", "2.5e2",
                              Patterns::Double(),
                              "Transport coefficient for shale in m^/2yr");
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
      end_time = prm.get_double ("End time");
      if (prm.get_bool ("Use years in output instead of seconds") == true)
        end_time *= year_in_seconds;
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection("Fastscape");
        {
          nstep = prm.get_integer("Number of steps");
          max_timestep = prm.get_double("Maximum timestep");
          vexp = prm.get_double("Vertical exaggeration");
          additional_refinement = prm.get_integer("Additional fastscape refinement");
          slice = prm.get_bool("Use center slice for 2d");
          fs_seed = prm.get_integer("Fastscape seed");
          surface_resolution = prm.get_integer("Surface resolution");
          resolution_difference = prm.get_integer("Resolution difference");
          use_marine = prm.get_bool("Use marine component");
          use_strat = prm.get_bool("Use stratigraphy");
          nstepp = prm.get_integer("Total steps");
          nreflectorp = prm.get_integer("Number of horizons");
          y_extent_2d = prm.get_double("Y extent in 2d");
          use_ghost = prm.get_bool("Use ghost nodes");
          use_v = prm.get_bool("Use velocities");
          precision = prm.get_double("Precision");
          sr_values = Utilities::string_to_double
                             (Utilities::split_string_list(prm.get ("Sediment rain")));
          sr_times = Utilities::string_to_double
                             (Utilities::split_string_list(prm.get ("Sediment rain intervals")));

              if (sr_values.size() != sr_times.size()+1)
                  AssertThrow(false, ExcMessage("Error: There must be one more sediment rain value than interval."));


          prm.enter_subsection("Boundary conditions");
          {
            bottom = prm.get_integer("Bottom");
            right = prm.get_integer("Right");
            top = prm.get_integer("Top");
            left = prm.get_integer("Left");
            left_flux = prm.get_double("Left mass flux");
            right_flux = prm.get_double("Right mass flux");
            top_flux = prm.get_double("Top mass flux");
            bottom_flux = prm.get_double("Bottom mass flux");

            // Put the boundary condition values into a four digit value to send to FastScape.
            bc = bottom*1000+right*100+top*10+left;

            if ((left_flux != 0 && top_flux != 0) || (left_flux != 0 && bottom_flux != 0) ||
                (right_flux != 0 && bottom_flux != 0) || (right_flux != 0 && top_flux != 0))
              AssertThrow(false,ExcMessage("Currently the plugin does not support mass flux through adjacent boundaries."));
          }
          prm.leave_subsection();

          prm.enter_subsection("Erosional parameters");
          {
            m = prm.get_double("Drainage area exponent");
            n = prm.get_double("Slope exponent");
            kfsed = prm.get_double("Sediment river incision rate");
            kff = prm.get_double("Bedrock river incision rate");
            kdsed = prm.get_double("Sediment diffusivity");
            kdd = prm.get_double("Bedrock diffusivity");
            g = prm.get_double("Bedrock deposition coefficient");
            gsed = prm.get_double("Sediment deposition coefficient");
            p = prm.get_double("Multi-direction slope exponent");
          }
          prm.leave_subsection();

          prm.enter_subsection("Marine parameters");
          {
            sl = prm.get_double("Sea level");
            p1 = prm.get_double("Sand porosity");
            p2 = prm.get_double("Shale porosity");
            z1 = prm.get_double("Sand e-folding depth");
            z2 = prm.get_double("Shale e-folding depth");
            r = prm.get_double("Sand-shale ratio");
            l = prm.get_double("Depth averaging thickness");
            kds1 = prm.get_double("Sand transport coefficient");
            kds2 = prm.get_double("Shale transport coefficient");
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
    ASPECT_REGISTER_MESH_DEFORMATION_MODEL(FastScape,
                                           "fastscape",
                                           "A plugin which uses the program FastScape to add surface processes "
                                           "such as diffusion, sedimentation, and the Stream Power Law to ASPECT. "
                                           "FastScape is initialized with ASPECT height and velocities, and then "
                                           "continues to run in the background, updating with new ASPECT velocities "
                                           "when called. Once FastScape has run for the given amount of timesteps, it "
                                           "compares the initial and final heights to send a velocity back to the mesh "
                                           "deformation plugin.")
  }
}
