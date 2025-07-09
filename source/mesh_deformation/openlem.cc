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



#include "aspect/global.h"
#include <cmath>
#include <deal.II/base/exceptions.h>
#include <deal.II/numerics/vector_tools.h>
#include <aspect/mesh_deformation/openlem.h>
#include <aspect/geometry_model/box.h>
#include <aspect/solution_evaluator.h>
#include <aspect/geometry_model/initial_topography_model/interface.h>
#include <limits>
#include <string>


namespace aspect
{
  namespace MeshDeformation
  {
    template <int dim>
    OpenLEM<dim>::OpenLEM()
      :
      connector(&grid_new)
    {}


    template <int dim>
    void
    OpenLEM<dim>::initialize ()
    {
      //if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
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
        const unsigned int aspect_refinement_level = maximum_surface_refinement_level;
        const unsigned int openlem_nodes = Utilities::pow(2,openlem_refinement_level);
        const unsigned int aspect_nodes = Utilities::pow(2,aspect_refinement_level);
        openlem_nx = openlem_nodes * repetitions[0] + ghost_nodes + 1;
        aspect_nx = aspect_nodes * repetitions[0] + ghost_nodes + 1;

        // Size of openlem cell.
        openlem_dx = (grid_extent[0].second)/(openlem_nodes * repetitions[0]);
        aspect_dx = (grid_extent[0].second)/(aspect_nodes * repetitions[0]);

        // openlem X extent, which is generally ASPECT's extent unless the ghost nodes are used,
        // in which case 2 cells are added on either side.
        openlem_x_extent = (grid_extent[0].second) + openlem_dx * ghost_nodes;
        aspect_x_extent = (grid_extent[0].second) + aspect_dx * ghost_nodes;

        // Sub intervals are 3 less than points, if including the ghost nodes. Otherwise 1 less.
        table_intervals[0] = openlem_nodes * repetitions[0];
        table_intervals[dim-1] = 1;

        if (dim == 2)
          {
            openlem_dy = openlem_dx;
            openlem_y_extent = std::round(openlem_y_extent_2d/openlem_dy)*openlem_dy + openlem_dy * ghost_nodes;
            openlem_ny = 1+openlem_y_extent/openlem_dy;
            aspect_dy = aspect_dx;
            aspect_y_extent = std::round(openlem_y_extent_2d/aspect_dy)*aspect_dy + aspect_dy * ghost_nodes;
            aspect_ny = 1+aspect_y_extent/aspect_dy;
          }
        else
          {
            openlem_ny = openlem_nodes * repetitions[1] + ghost_nodes + 1;
            openlem_dy = (grid_extent[1].second)/(openlem_nodes * repetitions[1]);
            table_intervals[1] = openlem_nodes * repetitions[1];
            openlem_y_extent = (grid_extent[1].second) + openlem_dy * ghost_nodes;
            aspect_ny = aspect_nodes * repetitions[1] + ghost_nodes + 1;
            aspect_dy = (grid_extent[1].second)/(aspect_nodes * repetitions[1]);
            openlem_y_extent = (grid_extent[1].second) + aspect_dy * ghost_nodes;
          }


        // Create a folder for the openlem visualization files.
        Utilities::create_directory (this->get_output_directory() + "openLEM/",
                                     this->get_mpi_communicator(),
                                     false);

        last_output_time = 0;
        // preparte the openLEM variables
        grid_old = openlem::OceanGrid(openlem_nx,openlem_ny);
        grid_new = openlem::OceanGrid(openlem_nx,openlem_ny);
        std::cout << "nx = " << openlem_nx << ", ny = " << openlem_ny << std::endl;

        connector = openlem::Connector(&grid_new,openlem_dy,grid_extent[0].first, grid_extent[1].first  );
        // todo: fill the grid new

        // Define all nodes with non-positive elevations (here, the ocean around the
        // island as boundary nodes
        for ( int i = 0; i < grid_new.m; ++i )
          for ( int j = 0; j < grid_new.n; ++j )
            {
              //grid_new[i][j].b = grid_new[i][j].h  <= 0 || i == 0 || j == 0 || i == grid_new.m-1 || j == grid_new.n-1;
              grid_new[i][j].b = i == 0 || j == 0 || i == grid_new.m-1 || j == grid_new.n-1
                                 || i == 1 || j == 1 || i == grid_new.m-2 || j == grid_new.n-2
                                 ;
              //grid_new[i][j].h = (grid_new[i][j].b =
              //                     i == 0 || j == 0 || i == grid_new.m-1 || j == grid_new.n-1) ? 1 : this->get_initial_topography_model().value(Point<dim-1>(connector.x[i][j],connector.y[i][j])); //||
              //i == 1 || j == 1 || i == grid_new.m-2 || j == grid_new.n-2
              //;
              grid_new[i][j].h = i == 0 || j == 0 || i == grid_new.m-1 || j == grid_new.n-1 || i == 1 || j == 1 || i == grid_new.m-2 || j == grid_new.n-2 ? 1 : this->get_initial_topography_model().value(Point<dim-1>(connector.x[i][j],connector.y[i][j]));
              //grid_new[i][j].b += (grid_new[i][j].h  <= 0)*2;
            }

        // Compute initial flow pattern and water level
        deepest_point = openlem::Point(0,0);
        for ( int i = 0; i < grid_new.m; ++i )
          for ( int j = 0; j < grid_new.n; ++j )
            {
              if ( grid_new[i][j].h < grid_new[deepest_point].h )
                deepest_point = openlem::Point(i,j);
            }
        grid_new[deepest_point].b = 1;
        grid_new.odiff = openlem_ocean_diffusivity;
        grid_new.fillLakes();
        grid_new.computeFluxes();
        grid_old = grid_new;
        std::cout << "opelem_dy = " << openlem_dy << ", dox = " << grid_extent[0].first << ", doy" << grid_extent[1].first << ", deepest point = " << deepest_point.i << ":" << deepest_point.j  << ", q of deepest point = " << grid_new[deepest_point].q << ", deepest point h = " << grid_new[deepest_point].h << ", b = " << (int)grid_new[deepest_point].b << std::endl;
        connector = openlem::Connector(&grid_new,openlem_dy,grid_extent[0].first, grid_extent[1].first, 0.5);
        //connector.interpolation(double &x, double &y)
        for (unsigned int x_i = 0; x_i < openlem_nx; ++x_i)
          for (unsigned int y_i= 0; y_i < openlem_ny; ++y_i)
            {
              //std::cout << "a x:y = " << x_i << ":" << y_i << " = " << connector.x[x_i][y_i] << ":" << connector.y[x_i][y_i] << std::endl;
              //connector.interpolation(connector.x[x_i][y_i], connector.y[x_i][y_i]);
            }
        connector = openlem::Connector(&grid_new,openlem_dy,grid_extent[0].first, grid_extent[1].first  );
        //for (unsigned int x_i = 0; x_i < openlem_nx; ++x_i)
        //  for (unsigned int y_i= 0; y_i < openlem_ny; ++y_i)
        //    std::cout << "x:y = " << x_i << ":" << y_i << " = " << connector.x[x_i][y_i] << ":" << connector.y[x_i][y_i] << std::endl;
      }
      openlem_mesh_aspect_cell = std::vector<std::vector<CellId>>(openlem_nx,std::vector<CellId>(openlem_ny));
    }


    template <int dim>
    void
    OpenLEM<dim>::update ()
    {
      std::string dirname = this->get_output_directory();
      const char *dirname_char=dirname.c_str();
      const unsigned int dirname_length = dirname.length();

      TimerOutput::Scope timer_section(this->get_computing_timer(), "openlem plugin");

      std::cout << "Flag 1: before (old:new)= " << (grid_old.getNode(5,7))->h << ":" << (grid_new.getNode(5,7))->h << std::endl;
      if (this->get_timestep_number() == 0)
        {
          for (unsigned int x_i=0; x_i<openlem_nx; ++x_i)
            for (unsigned int y_i=0; y_i<openlem_ny; ++y_i)
              {
                //openlem::Node *node = grid_new.getNode(x_i,y_i);
                // TODO: should I use the openlem location here or the aspect projected point?
                //node->h = this->get_initial_topography_model().value(Point<dim-1>(connector.x[x_i][y_i],connector.y[x_i][y_i]));
              }
          if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            connector.write_vtk(grid_extent[dim-1].second, 0, 0., dirname);
          return;
        }
      grid_old = grid_new;
      //connector = openlem::Connector(&grid_new,openlem_dy);
      std::cout << "Flag 2 before (old:new)= " << (grid_old.getNode(5,7))->h << ":" << (grid_new.getNode(5,7))->h << std::endl;

      const unsigned int current_timestep = this->get_timestep_number ();
      const double aspect_timestep_in_years = this->get_timestep() / year_in_seconds;

      // Find a openlem timestep that is below our maximum timestep.
      unsigned int openlem_iterations = openlem_steps_per_aspect_step;
      double openlem_timestep_in_years = aspect_timestep_in_years/openlem_iterations;
      //while (openlem_timestep_in_years>maximum_openlem_timestep)
      {
        //openlem_iterations *= 2;
        //openlem_timestep_in_years *= 0.5;
      }

      // Vector to hold the velocities that represent the change to the surface.
      //const unsigned int openlem_array_size = openlem_nx*openlem_ny;
      std::vector<std::vector<Point<dim>>> mesh_locations(openlem_nx,std::vector<Point<dim>>(openlem_ny));
      // fill mesh locations
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          for (unsigned int x_i = 0; x_i < openlem_nx; ++x_i)
            {
              for (unsigned int y_i = 0; y_i < openlem_nx; ++y_i)
                {
                  // The z will be replaced later on
                  if (dim == 3)
                    {
                      mesh_locations[x_i][y_i] = Point<dim>(connector.x[x_i][y_i],connector.y[x_i][y_i],grid_extent[dim-1].second + grid_new[x_i][y_i].h);
                      //if (x_i == 10)
                      //  std::cout << "A: rank = " <<Utilities::MPI::this_mpi_process(this->get_mpi_communicator())  << ", x:y = " << x_i << ":" << y_i << ", mesh_locations[x_i][y_i] = " << mesh_locations[x_i][y_i] << ", connector x:y = " <<  connector.x[x_i][y_i] << ":" << connector.x[x_i][y_i] << std::endl;
                      // NOTE: at timestep 0
                      //if (this->get_timestep_number() == 0)
                      //  {
                      //    // TODO: probably should move this if statement out of the double loop
                      //    //std::cout << "it: " << this->get_initial_topography_model().value(Point<dim-1>(connector.x[x_i][y_i],connector.y[x_i][y_i])) << ", x = " << connector.x[x_i][y_i] << ", y = " << connector.y[x_i][y_i] << std::endl;
                      //    mesh_locations[x_i][y_i] = Point<dim>(connector.x[x_i][y_i],connector.y[x_i][y_i],grid_extent[dim-1].second);// + this->get_initial_topography_model().value(Point<dim-1>(connector.x[x_i][y_i],connector.y[x_i][y_i])));
                      //  }
                      //else
                      //  {
                      //    mesh_locations[x_i][y_i] = Point<dim>(connector.x[x_i][y_i],connector.y[x_i][y_i],grid_extent[dim-1].second);// + grid_new[x_i][y_i].h);
                      //  }
                    }
                }
            }
        }


      std::cout << "before  broadcast mesh_locations" << std::endl;
      mesh_locations = Utilities::MPI::broadcast(this->get_mpi_communicator(),mesh_locations,0);
      std::cout << "after  broadcast mesh_locations" << std::endl;

      for (unsigned int x_i = 0; x_i < openlem_nx; ++x_i)
        {
          for (unsigned int y_i = 0; y_i < openlem_nx; ++y_i)
            {
              if (dim == 3)
                {
                  //if (x_i == 10)
                  //  std::cout << "B: rank = " <<Utilities::MPI::this_mpi_process(this->get_mpi_communicator())  << ", x:y = " << x_i << ":" << y_i << ", mesh_locations[x_i][y_i] = " << mesh_locations[x_i][y_i] << ", connector x:y = " <<  connector.x[x_i][y_i] << ":" << connector.y[x_i][y_i] << std::endl;
                }
            }
        }
      std::vector<std::vector<Tensor<1,dim>>> mesh_velocities(openlem_nx,std::vector<Tensor<1,dim>>(openlem_ny));
      std::vector<std::vector<Point<dim>>> mesh_velocity_locations(openlem_nx,std::vector<Point<dim>>(openlem_ny));
      std::vector<std::vector<double>> mesh_velocity_distances(openlem_nx,std::vector<double>(openlem_ny));
      aspect_mesh_z = std::vector<std::vector<double>>(aspect_nx,std::vector<double>(aspect_ny));


      for (unsigned int x_i = 0; x_i < aspect_nx; ++x_i)
        {
          for (unsigned int y_i = 0; y_i < aspect_ny; ++y_i)
            {
              aspect_mesh_z[x_i][y_i] = std::numeric_limits<double>::signaling_NaN();
            }
        }
      for (unsigned int x_i = 0; x_i < openlem_nx; ++x_i)
        {
          for (unsigned int y_i = 0; y_i < openlem_ny; ++y_i)
            {
              mesh_velocities[x_i][y_i][0] = std::numeric_limits<double>::signaling_NaN();
              mesh_velocities[x_i][y_i][1] = std::numeric_limits<double>::signaling_NaN();
              mesh_velocities[x_i][y_i][dim-1] = std::numeric_limits<double>::signaling_NaN();
              mesh_velocity_locations[x_i][y_i][0] = std::numeric_limits<double>::signaling_NaN();
              mesh_velocity_locations[x_i][y_i][1] = std::numeric_limits<double>::signaling_NaN();
              mesh_velocity_locations[x_i][y_i][dim-1] = std::numeric_limits<double>::signaling_NaN();
              mesh_velocity_distances[x_i][y_i] = std::numeric_limits<double>::infinity();
            }
        }
      std::vector<std::tuple<unsigned int, unsigned int, double,Point<dim>, Tensor<1,dim>>> mesh_velocities_vector;

      //const UpdateFlags update_flags = update_values | update_gradients;
      //std::unique_ptr<SolutionEvaluator<dim>> evaluator = construct_solution_evaluator(*this,update_flags);
      //std::vector<EvaluationFlags::EvaluationFlags> evaluation_flags (this->introspection().n_components, EvaluationFlags::nothing);

      //for (unsigned int i=0; i<this->introspection().n_components; ++i)
      //  {
      //    evaluation_flags[i] |= EvaluationFlags::values;
      //    evaluation_flags[i] |= EvaluationFlags::gradients;
      //  }



      const UpdateFlags update_flags = update_values;// | update_gradients;
      std::unique_ptr<SolutionEvaluator<dim>> evaluator = construct_solution_evaluator(*this,update_flags);

      // Todo: assert that vectors of cells, positions and reference positino are the same size

      // FEPointEvaluation uses different evaluation flags than the common UpdateFlags.
      // Translate between the two.
      std::vector<EvaluationFlags::EvaluationFlags> evaluation_flags (this->introspection().n_components, EvaluationFlags::nothing);

      for (unsigned int i=0; i<this->introspection().n_components; ++i)
        {
          evaluation_flags[i] |= EvaluationFlags::values;
          //evaluation_flags[i] |= EvaluationFlags::gradients;
        }

      std::vector<Tensor< 1, dim, double >> velocity_field;
      std::vector<std::vector<double>> solution(this->get_fe().dofs_per_cell);

      //std::vector<std::vector<Tensor<1,dim>>> gradients(this->get_fe().dofs_per_cell);
      solution.resize(1,std::vector<double>(evaluator->n_components(), numbers::signaling_nan<double>()));
      //gradients.resize(1,std::vector<Tensor<1,dim>>(evaluator->n_components(), numbers::signaling_nan<Tensor<1,dim>>()));

      const types::boundary_id relevant_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id ("top");
      // Get a quadrature rule that exists only on the corners, and increase the refinement if specified.
      const QIterated<dim-1> face_corners (QTrapezoid<1>(),
                                           1);//additional_refinement_levels+surface_refinement_difference));

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        face_corners,
                                        update_values |
                                        update_quadrature_points);
      // first check if the stored cell is still the cell we need.
      //std::cout << "face_corners.size() = " << face_corners.size() << std::endl;
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned() && cell->at_boundary())
          for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
            if (cell->face(face_no)->at_boundary())
              {
                if ( cell->face(face_no)->boundary_id() != relevant_boundary)
                  continue;

                //std::vector<Tensor<1,dim>> vel(face_corners.size());
                fe_face_values.reinit(cell, face_no);
                //fe_face_values[this->introspection().extractors.velocities].get_function_values(this->get_solution(), vel);

                for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
                  {
                    const Point<dim> vertex = fe_face_values.quadrature_point(corner);

                    const double x_coord = vertex[0] - grid_extent[0].first;
                    const double y_coord = vertex[1] - grid_extent[1].first;
                    const size_t ixa = std::round(x_coord/aspect_dx);
                    const size_t iya = std::round(y_coord/aspect_dy);
                    aspect_mesh_z[ixa][iya] = vertex[dim-1] - grid_extent[dim-1].second;
                  }
              }

      //unsigned int total_found_velocities = 0;
      const CellId unitialized_cell_id = CellId();
      for (unsigned int x_i = 0; x_i < openlem_nx; ++x_i)
        {
          for (unsigned int y_i = 0; y_i < openlem_ny; ++y_i)
            {
              bool found_a_value = false;
              std::vector<Point<dim>> unit_point(1);
              double closest_distance_cell = std::numeric_limits<double>::infinity();
              //std::cout << "cellid = " << openlem_mesh_aspect_cell[x_i][y_i].get_coarse_cell_id() << "..." << std::endl;
              //std::cout << "try cell id at " << x_i << ":" << y_i << ", coarse cell id = " <<  openlem_mesh_aspect_cell[x_i][y_i].get_coarse_cell_id()<< std::endl;
              if (openlem_mesh_aspect_cell[x_i][y_i].get_coarse_cell_id() != numbers::invalid_coarse_cell_id)
                {
                  //std::cout << "found a valid cell id at " << x_i << ":" << y_i << std::endl;
                  const auto &cell = this->get_triangulation().create_cell_iterator(openlem_mesh_aspect_cell[x_i][y_i])->as_dof_handler_iterator(this->get_dof_handler());
                  if ( cell->state() == IteratorState::valid && cell->is_locally_owned() && cell->at_boundary())
                    {
                      //std::cout << "found a valid cell at " << x_i << ":" << y_i << std::endl;
                      for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
                        {
                          if (cell->face(face_no)->at_boundary())
                            {
                              if ( cell->face(face_no)->boundary_id() != relevant_boundary)
                                continue;

                              this->get_mapping().transform_points_real_to_unit_cell(cell, {mesh_locations[x_i][y_i]},unit_point);
                              const double distance_to_unit_cell = GeometryInfo<dim>::distance_to_unit_cell(unit_point[0]);

                              if (distance_to_unit_cell < closest_distance_cell)
                                {
                                  std::vector<Point<dim>> closest_point = {cell->reference_cell().closest_point(unit_point[0])};
                                  Point<dim> closest_point_real = this->get_mapping().transform_unit_to_real_cell(cell,closest_point[0]);

                                  Point<2> cp2d(closest_point_real[0],closest_point_real[1]);
                                  Point<2> openlem_point_2d(connector.x[x_i][y_i],connector.y[x_i][y_i]);


                                  if (cell->state() != IteratorState::valid)
                                    {
                                      continue;
                                    }

                                  auto &velocity_evaluator = evaluator->get_velocity_or_fluid_velocity_evaluator(false);
                                  small_vector<double,50> solution_values(this->get_fe().dofs_per_cell);
                                  cell->get_dof_values(this->get_current_linearization_point(),//input_solution,
                                                       solution_values.begin(),
                                                       solution_values.end());


                                  solution[0] = std::vector<double>(evaluator->n_components(), numbers::signaling_nan<double>());


                                  evaluator->reinit(cell, closest_point);
                                  //evaluator->evaluate({solution_values.data(),solution_values.size()},evaluation_flags);
                                  //evaluator->get_solution(0, {&solution[0][0],solution[0].size()}, evaluation_flags);

                                  velocity_evaluator.evaluate({solution_values.data(),solution_values.size()},
                                                              EvaluationFlags::values);

                                  auto &mapping_info = evaluator->get_mapping_info();
                                  mapping_info.reinit(cell, {closest_point.data(),closest_point.size()});


                                  Point<3> velocity_3d;
                                  Point<3> velocity_original;
                                  Point<3> velocity_analytical(-openlem_point_2d[1],openlem_point_2d[0],0.0);
                                  Tensor<1,dim> velocity = velocity_evaluator.get_value(0)*year_in_seconds;

                                  //if (y_i == 10)
                                  //  std::cout << "rank = " <<Utilities::MPI::this_mpi_process(this->get_mpi_communicator())  << ", x:y = " << x_i << ":" << y_i << ", closest_point = " << closest_point[0] << ", closest_point_real = " << closest_point_real<< ", mesh_locations = " << mesh_locations[x_i][y_i] << ", velo = " <<  velocity << ", velo ana = " << velocity_analytical << ", distance_to_unit_cell = " << distance_to_unit_cell<< ", vel raw = " << velocity_evaluator.get_value(0) << std::endl;
                                  for (unsigned int i = 0; i < dim; ++i)
                                    {
                                      velocity_original[i] = velocity[i];
                                      //velocity[i] = velocity_analytical[i];// * year_in_seconds;//solution_values[this->introspection().component_indices.velocities[i]] * year_in_seconds;
                                      velocity_3d[i] = velocity[i];//solution_values[this->introspection().component_indices.velocities[i]] * year_in_seconds;
                                    }

                                  // TODO: maybe don't need the closest point real in the tuple?
                                  //std::tuple<unsigned int, unsigned int, Point<dim>, Tensor<1,dim>> tuple = std::make_tuple(x_i, y_i, closest_point_real, velocity);
                                  if (!isnan(velocity[0]))
                                    {
                                      found_a_value = true;
                                      closest_distance_cell = distance_to_unit_cell; // TODO: probably just skip out if distance is exactly zero.
                                      //total_found_velocities++;
                                      mesh_velocities_vector.emplace_back(std::make_tuple(x_i, y_i, distance_to_unit_cell, closest_point_real, velocity));
                                      openlem_mesh_aspect_cell[x_i][y_i] = cell->id();
                                      //std::cout << "set cell id at " << x_i << ":" << y_i << ", coarse cell id = " <<  cell->id().get_coarse_cell_id() << ", " << openlem_mesh_aspect_cell[x_i][y_i]  << std::endl;
                                      break;
                                    }
                                }
                            }
                        }
                    }
                }
              if (closest_distance_cell != 0.0)
                {
                  for (const auto &cell : this->get_dof_handler().active_cell_iterators())
                    {
                      if (closest_distance_cell == 0.0)
                        {
                          // if the point is inside the cell it is defined as exactly zero
                          break;
                        }
                      if (cell->is_locally_owned() && cell->at_boundary())
                        {
                          for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
                            {
                              if (cell->face(face_no)->at_boundary())
                                {
                                  if ( cell->face(face_no)->boundary_id() != relevant_boundary)
                                    continue;

                                  this->get_mapping().transform_points_real_to_unit_cell(cell, {mesh_locations[x_i][y_i]},unit_point);
                                  const double distance_to_unit_cell = GeometryInfo<dim>::distance_to_unit_cell(unit_point[0]);


                                  if (distance_to_unit_cell < closest_distance_cell)
                                    {
                                      std::vector<Point<dim>> closest_point = {cell->reference_cell().closest_point(unit_point[0])};
                                      Point<dim> closest_point_real = this->get_mapping().transform_unit_to_real_cell(cell,closest_point[0]);

                                      Point<2> cp2d(closest_point_real[0],closest_point_real[1]);
                                      Point<2> openlem_point_2d(connector.x[x_i][y_i],connector.y[x_i][y_i]);


                                      if (cell->state() != IteratorState::valid)
                                        {
                                          continue;
                                        }

                                      auto &velocity_evaluator = evaluator->get_velocity_or_fluid_velocity_evaluator(false);
                                      small_vector<double,50> solution_values(this->get_fe().dofs_per_cell);
                                      cell->get_dof_values(this->get_current_linearization_point(),//input_solution,
                                                           solution_values.begin(),
                                                           solution_values.end());


                                      solution[0] = std::vector<double>(evaluator->n_components(), numbers::signaling_nan<double>());


                                      evaluator->reinit(cell, closest_point);
                                      //evaluator->evaluate({solution_values.data(),solution_values.size()},evaluation_flags);
                                      //evaluator->get_solution(0, {&solution[0][0],solution[0].size()}, evaluation_flags);

                                      velocity_evaluator.evaluate({solution_values.data(),solution_values.size()},
                                                                  EvaluationFlags::values);

                                      auto &mapping_info = evaluator->get_mapping_info();
                                      mapping_info.reinit(cell, {closest_point.data(),closest_point.size()});


                                      Point<3> velocity_3d;
                                      Point<3> velocity_original;
                                      Point<3> velocity_analytical(-openlem_point_2d[1],openlem_point_2d[0],0.0);
                                      Tensor<1,dim> velocity = velocity_evaluator.get_value(0)*year_in_seconds;

                                      //if (y_i == 10)
                                      //  std::cout << "rank = " <<Utilities::MPI::this_mpi_process(this->get_mpi_communicator())  << ", x:y = " << x_i << ":" << y_i << ", closest_point = " << closest_point[0] << ", closest_point_real = " << closest_point_real<< ", mesh_locations = " << mesh_locations[x_i][y_i] << ", velo = " <<  velocity << ", velo ana = " << velocity_analytical << ", distance_to_unit_cell = " << distance_to_unit_cell<< ", vel raw = " << velocity_evaluator.get_value(0) << std::endl;
                                      for (unsigned int i = 0; i < dim; ++i)
                                        {
                                          velocity_original[i] = velocity[i];
                                          //velocity[i] = velocity_analytical[i];// * year_in_seconds;//solution_values[this->introspection().component_indices.velocities[i]] * year_in_seconds;
                                          velocity_3d[i] = velocity[i];//solution_values[this->introspection().component_indices.velocities[i]] * year_in_seconds;
                                        }

                                      // TODO: maybe don't need the closest point real in the tuple?
                                      //std::tuple<unsigned int, unsigned int, Point<dim>, Tensor<1,dim>> tuple = std::make_tuple(x_i, y_i, closest_point_real, velocity);
                                      if (!isnan(velocity[0]))
                                        {
                                          found_a_value = true;
                                          closest_distance_cell = distance_to_unit_cell; // TODO: probably just skip out if distance is exactly zero.
                                          //total_found_velocities++;
                                          mesh_velocities_vector.emplace_back(std::make_tuple(x_i, y_i, distance_to_unit_cell, closest_point_real, velocity));
                                          openlem_mesh_aspect_cell[x_i][y_i] = cell->id();
                                          //std::cout << "set cell id at " << x_i << ":" << y_i << ", coarse cell id = " <<  cell->id().get_coarse_cell_id() << ", " << openlem_mesh_aspect_cell[x_i][y_i]  << std::endl;
                                          break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
              else
                {
                  //std::cout << "Found the right cell in storage! " << std::endl;
                }
              AssertThrow(found_a_value,ExcMessage("Could not found a value for " + std::to_string(x_i) + ":" + std::to_string(y_i)));
            }
        }

      // the mesh velocties vector should contain all the entries now, return to process 0 and unpack
      std::cout << "before gather mesh_velocites_vector" << std::endl;
      std::vector<std::vector<std::vector<double>>> aspect_mesh_z_vectors = Utilities::MPI::gather(this->get_mpi_communicator(),aspect_mesh_z,0);
      std::vector<std::vector<std::tuple<unsigned int, unsigned int, double,Point<dim>, Tensor<1,dim>>>> mesh_velocties_vectors = Utilities::MPI::gather(this->get_mpi_communicator(),mesh_velocities_vector,0);
      std::cout << "after gather mesh_velocites_vector" << std::endl;
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          for (auto &aspect_mesh_z_item : aspect_mesh_z_vectors)
            {
              for (size_t xi = 0; xi < aspect_nx; ++xi)
                {
                  for (size_t yi = 0; yi < aspect_ny; ++yi)
                    {

                      if (!std::isnan(aspect_mesh_z_item[xi][yi]))
                        {
                          aspect_mesh_z[xi][yi] = aspect_mesh_z_item[xi][yi];
                        }
                    }

                }
            }
          size_t counter = 0;
          for (auto &mesh_velocities_vector : mesh_velocties_vectors)
            {
              for (auto &position_velocity : mesh_velocities_vector)
                {
                  //if (std::get<0>(position_velocity) == 10)
                  //  std::cout << counter << "a: rank = " <<Utilities::MPI::this_mpi_process(this->get_mpi_communicator())  << ", mesh velocties x:y = " << std::get<0>(position_velocity) << ":" << std::get<1>(position_velocity) << ", velo = " << std::get<3>(position_velocity) << ", distance = " <<  std::get<2>(position_velocity)<< ", mesh_velocity_distances  = " << mesh_velocity_distances[std::get<0>(position_velocity)][std::get<1>(position_velocity)] << std::endl;
                  if (!std::isnan(std::get<4>(position_velocity)[0]) && std::get<2>(position_velocity) < mesh_velocity_distances[std::get<0>(position_velocity)][std::get<1>(position_velocity)] )
                    {
                      //   if (std::get<0>(position_velocity) == 10)
                      //     std::cout << counter << "b: rank = " <<Utilities::MPI::this_mpi_process(this->get_mpi_communicator())  << ", mesh velocties x:y = " << std::get<0>(position_velocity) << ":" << std::get<1>(position_velocity) << ", velo = " << std::get<3>(position_velocity) << ", distance = " <<  std::get<2>(position_velocity)<< std::endl;
                      mesh_velocities[std::get<0>(position_velocity)][std::get<1>(position_velocity)] = std::get<4>(position_velocity);
                      mesh_velocity_locations[std::get<0>(position_velocity)][std::get<1>(position_velocity)] = std::get<3>(position_velocity);
                      mesh_velocity_distances[std::get<0>(position_velocity)][std::get<1>(position_velocity)] = std::get<2>(position_velocity);
                    }
                }
              counter++;
            }
          // TODO: check whether the x and y of the closest_point_real and the original connector point are actually the same (or close enought)


          //std::vector<std::vector<double>> local_aspect_values = get_aspect_values();
          // Because there is no increase in time during timestep 0, we return and only
          // initialize and run openlem from timestep 1 and on.
          for (unsigned int x_i=0; x_i<openlem_nx; ++x_i)
            for (unsigned int y_i=0; y_i<openlem_ny; ++y_i)
              {
                openlem::Node *node = grid_new.getNode(x_i,y_i);
                // TODO: should I use the openlem location here or the aspect projected point?
                if (this->get_timestep_number() == 0)
                  {
                    node->h = this->get_initial_topography_model().value(Point<dim-1>(connector.x[x_i][y_i],connector.y[x_i][y_i]));
                  }

                node->l = 0;
                node->u = mesh_velocities[x_i][y_i][dim-1];//*100000;//1;//local_aspect_values[dim+2][i];
                // TODO: 2D
                connector.vx[x_i][y_i] =  mesh_velocities[x_i][y_i][0];//*100000;
                connector.vy[x_i][y_i] =  mesh_velocities[x_i][y_i][1];//*100000;
                connector.vz[x_i][y_i] =  mesh_velocities[x_i][y_i][2];//*100000;
              }
          //connector.write_vtk(grid_extent[dim-1].second,this->get_timestep_number(), 0.0, dirname  , "_preconvert");
          connector.convertVelocities(openlem_minimize_advection);
          //connector.write_vtk(grid_extent[dim-1].second, this->get_timestep_number() , "postconvert");
          connector.computeUpliftRate();
          //connector.write_vtk(grid_extent[dim-1].second, this->get_timestep_number() , "postuplift");
          // Define all nodes with non-positive elevations (here, the ocean around the
          // island as boundary nodes
          for ( int i = 0; i < grid_new.m; ++i )
            for ( int j = 0; j < grid_new.n; ++j )
              {
                //grid_new[i][j].b = grid_new[i][j].h <= 0;
                //  grid_new[i][j].b =  i == 0 || j == 0 || i == grid_new.m-1 || j == grid_new.n;
                //grid_new[i][j].b = 0;
                //i == 0 || j == 0 || i == grid_new.m-1 || j == grid_new.n-1 ||
                //i == 1 || j == 1 || i == grid_new.m-2 || j == grid_new.n-2
                //;
                //grid_new[i][j].b += (grid_new[i][j].h  <= 0)*2;
                //    grid_new[i][j].u = grid_new[i][j].h <= 0;
              }
          openlem::OceanGrid &g = grid_new;
          grid_new.kt = openlem_kt;//1e-14;
          grid_new.kd = openlem_kd;// 1e-14;

          //grid_new.fillLakes();
          // find deepest point
          // call set deepest point to 1
          // fill lakes
          // in each timestep compute water level()
          // call function clear ocean()
          // call function mark ocean ()
          //

        }
      if (this->get_timestep_number() == 0)
        {
          connector.write_vtk(grid_extent[dim-1].second, 0, 0.0, dirname);
          return;
        }


      // Run openlem on single process.
      //mesh_velocity_z = std::vector<std::vector<double>>(openlem_nx,std::vector<double>(openlem_ny));
      aspect_mesh_dh = std::vector<std::vector<double>>(aspect_nx,std::vector<double>(aspect_ny));
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          // Initialize the variables that will be sent to openlem.
          // Elevation is initialized at a very high number so that we can later check that all points
          // received data from ASPECT, and if not throw an assert.

          //std::cout << "Flag 3 before (old:new)= " << (grid_old.getNode(5,7))->h << ":" << (grid_new.getNode(5,7))->h << std::endl;

          //for (unsigned int ix=0; ix<openlem_nx; ++ix)
          //  for (unsigned int iy=0; iy<openlem_ny; ++iy)
          //    {
          //      openlem::Node *node = grid.getNode(index_x,index_y);
          //      node->h = local_aspect_values[0][i];
          //      node->l = 0;
          //      node->u = ;//1;//local_aspect_values[dim+2][i];
          //    }
          //fill_openlem_arrays(grid_new,
          //                    local_aspect_values);
          //std::cout << "Flag 4 before (old:new)= " << (grid_old.getNode(5,7))->h << ":" << (grid_new.getNode(5,7))->h << std::endl;

          if (current_timestep == 1 || restart)
            {
              this->get_pcout() << "   Initializing openlem... " << (1+maximum_surface_refinement_level+additional_refinement_levels) <<
                                " levels, cell size: " << openlem_dx << " m." << std::endl;

              // TODO
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
          const double time_in_years = this->get_time();// / year_in_seconds;
          //auto it = std::lower_bound(sediment_rain_times.begin(), sediment_rain_times.end(), time_in_years);
          //const unsigned int inds = std::distance(sediment_rain_times.begin(), it);
          //const double sediment_rain = sediment_rain_rates[inds];

          grid_old = grid_new;
          std::cout << "main before (old:new)= " << (grid_old.getNode(5,7))->h << ":" << (grid_new.getNode(5,7))->h << std::endl;
          execute_openlem(grid_new,openlem_timestep_in_years,openlem_iterations,dirname);

          //connector.updateCoordinateSystem (aspect_timestep_in_years);
          std::cout << "main after (old:new)= " << (grid_old.getNode(5,7))->h << ":" << (grid_new.getNode(5,7))->h << std::endl;
          // Write a file to store h, b & step for restarting.
          // TODO: It would be good to roll this into the general ASPECT checkpointing,
          // and when we do this needs to be changed.
          //if (((this->get_parameters().checkpoint_time_secs == 0) &&
          //     (this->get_parameters().checkpoint_steps > 0) &&
          //     ((current_timestep + 1) % this->get_parameters().checkpoint_steps == 0)) ||
          //    (this->get_time() == this->get_end_time() && this->get_timestepping_manager().need_checkpoint_on_terminate()))
          //  {
          //    //save_restart_files(elevation,
          //    //                   basement,
          //    //                   silt_fraction);
          //  }

          // Find out our velocities from the change in height.
          // Where mesh_velocity_z is a vector of array size that exists on all processes.
          //for (unsigned int ix=0; ix<openlem_nx; ++ix)
          //  for (unsigned int iy=0; iy<openlem_ny; ++iy)
          //    {
          //      openlem::Node *node_old = grid_old.getNode(ix,iy);
          //      openlem::Node *node_new = grid_new.getNode(ix,iy);
          //      mesh_velocity_z[ix][iy] = (node_new->h - node_old->h)/(1.0*this->get_timestep());//aspect_timestep_in_years;
          //      if ((ix == 18 || ix == 19) && iy == 29)
          //        std::cout << ix << ":" << iy << " = " << mesh_velocity_z[ix][iy] << ", node_old->h = " << node_old->h << ", node_new->h = " << node_new->h << ", node_new->h - node_old->h = " << node_new->h - node_old->h  << std::endl;
          //    }

          //Utilities::MPI::broadcast(this->get_mpi_communicator(), mesh_velocity_z, 0);


          // loop over all openlem grid points, find the closest aspect grid point
          // add the velocity to the closest grid point and count how many grid points contributed to that point
          // then divide the value by the counter.
          std::vector<std::vector<size_t>> aspect_mesh_counter = std::vector<std::vector<size_t>>(aspect_nx,std::vector<size_t>(aspect_ny));
          //double min_openlem_h = std::numeric_limits<double>::infinity();
          //double max_openlem_h = -std::numeric_limits<double>::infinity();
          for (unsigned int ix=0; ix<openlem_nx; ++ix)
            for (unsigned int iy=0; iy<openlem_ny; ++iy)
              {
                double x_coord = connector.x[ix][iy] - grid_extent[0].first;
                double y_coord = connector.y[ix][iy] - grid_extent[1].first;
                if (
                  x_coord >= 0 && x_coord <= grid_extent[0].second &&
                  y_coord >= 0 && y_coord <= grid_extent[1].second
                  //connector.x[ix][iy] >= grid_extent[0].first && connector.x[ix][iy] <= grid_extent[0].first+ grid_extent[0].second &&
                  //connector.y[ix][iy] >= grid_extent[1].first && connector.x[ix][iy] <= grid_extent[1].first+ grid_extent[1].second
                )
                  {

                    //Point<2> openlem_point (connector.x[ix][iy],connector.y[ix][iy]);
                    //Point<2> aspect_point_round (std::round((connector.x[ix][iy]-grid_extent[0].first)/aspect_dx),std::round((connector.y[ix][iy]-grid_extent[1].first)/aspect_dy));
                    // TODO: add assert that points are > 0
                    size_t closest_point_ixa = std::round(x_coord/aspect_dx);//aspect_point_round[0];//std::numeric_limits<size_t>::signaling_NaN();
                    size_t closest_point_iya = std::round(y_coord/aspect_dy);//aspect_point_round[1];//std::numeric_limits<size_t>::signaling_NaN();
                    //if (grid_new.getNode(ix,iy)->b == 0)
                    //min_openlem_h = std::min(min_openlem_h,grid_new.getNode(ix,iy)->h);
                    //max_openlem_h = std::max(max_openlem_h,grid_new.getNode(ix,iy)->h);
                    aspect_mesh_dh[closest_point_ixa][closest_point_iya] += grid_new.getNode(ix,iy)->h - aspect_mesh_z[closest_point_ixa][closest_point_iya];// + grid_extent[dim-1].second;//- mesh_velocity_locations[ix][iy][dim-1]+grid_extent[dim-1].second ;// grid_old.getNode(ix,iy)->h;//TODO: check if this is alright or if I need to get the closest node in the old grid sepeartly. //mesh_velocity_z[ix][iy];
                    //if (grid_new.getNode(ix,iy)->h > 10 )
                    //  std::cout << "rank = " <<  Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) << ", ix:iy = " << ix << ":" << iy << ", ixa:iya = " << closest_point_ixa << ":" << closest_point_iya << ", h = " << grid_new.getNode(ix,iy)->h << ", aspect_mesh_z =" << aspect_mesh_z[closest_point_ixa][closest_point_iya] << ", result = " << aspect_mesh_dh[closest_point_ixa][closest_point_iya] << std::endl;
                    aspect_mesh_counter[closest_point_ixa][closest_point_iya] += 1;
                  }
              }

          for (unsigned int ixa=0; ixa<aspect_nx; ++ixa)
            for (unsigned int iya=0; iya<aspect_ny; ++iya)
              {
                if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
                  {
                    //std::cout << "ixa:iya = " << ixa << ":" << iya << ",aspect_mesh_elevation_h = " << aspect_mesh_elevation_h[ixa][iya] << ", counter = " <<  aspect_mesh_counter[ixa][iya]<< std::endl;
                  }
                //if ( grid_new.getNode(ix,iy)->b == 0)
                {
                  if (aspect_mesh_counter[ixa][iya] == 0)
                    {
                      aspect_mesh_dh[ixa][iya] = 0;
                    }
                  else
                    {
                      AssertThrow(aspect_mesh_counter[ixa][iya] != 0, ExcMessage("Division by zero when averaging aspect mesh velocity z"));
                      //const double before = aspect_mesh_elevation_h[ixa][iya];
                      aspect_mesh_dh[ixa][iya] /= (aspect_mesh_counter[ixa][iya]);
                      //AssertThrow(aspect_mesh_elevation_h[ixa][iya] <= max_openlem_h && aspect_mesh_elevation_h[ixa][iya] >= min_openlem_h  , ExcMessage("wrong!!!! min:max h = " + std::to_string(min_openlem_h) + ":" + std::to_string(max_openlem_h) + ", value:after:counter = " + std::to_string(aspect_mesh_elevation_h[ixa][iya]) + ":" + std::to_string(before)   + ":" +  std::to_string(aspect_mesh_counter[ixa][iya])));
                    }
                  //if ((ixa == 18*2 || ixa == 19*2) && iya == 29*2)
                  //std::cout << "ix:iy = " << ixa << ":" << iya << ", aspect_mesh_elevation_h = " <<  aspect_mesh_elevation_h[ixa][iya] << std::endl;
                }

              }
        }
      else
        {
          std::cout << "main before (old:new) " << std::endl;
          execute_openlem(grid_new,openlem_timestep_in_years,openlem_iterations, dirname);

          std::cout << "main after (old:new) " << std::endl;
          //for (unsigned int i=0; i<local_aspect_values.size(); ++i)
          //MPI_Ssend(&local_aspect_values[i][0], local_aspect_values[1].size(), MPI_DOUBLE, 0, 42, this->get_mpi_communicator());

          // Check whether the openlem mesh was filled with data.
          // const bool openlem_mesh_filled = Utilities::MPI::broadcast (this->get_mpi_communicator(), true, 0);
          // if (openlem_mesh_filled != true)
          //   throw aspect::QuietException();

          // // This is called solely so we can set the timer and will return immediately.
          // execute_openlem(grid_new,
          //                 //mesh_velocity_z,
          //                 //mesh_velocity_z,
          //                 //mesh_velocity_z,
          //                 //mesh_velocity_z,
          //                 //mesh_velocity_z,
          //                 aspect_timestep_in_years,
          //                 openlem_steps_per_aspect_step);

          //mesh_velocity_z = Utilities::MPI::broadcast(this->get_mpi_communicator(), mesh_velocity_z, 0);
        }
      //std::cout << "before  broadcast mesh_velocity_z" << std::endl;
      aspect_mesh_dh = Utilities::MPI::broadcast(this->get_mpi_communicator(), aspect_mesh_dh, 0);
      //std::cout << "after  broadcast mesh_velocity_z, aspect_dx:dy = " << aspect_dx << ":" << aspect_dy << ", ge = " << grid_extent[0].first << ":" << grid_extent[1].first << ", openlem_nx = " << openlem_nx << ", aspect_nx = " << aspect_nx << std::endl;
      //std::cout << "min:max h = " << min_openlem_h << ":" << max_openlem_h << std::endl;
      // Get the sizes needed for a data table of the mesh velocities.
      //TableIndices<dim> size_idx;
      //for (unsigned int d=0; d<dim; ++d)
      //  {
      //    size_idx[d] = table_intervals[d]+1;
      //  }

      //std::cout << "table size = " << size_idx[0] << ":" << size_idx[1] << ":" << size_idx[2] << std::endl;
      // Initialize a table to hold all velocity values that will be interpolated back to ASPECT.
      //const Table<dim,double> velocity_table;// = fill_data_table(mesh_velocity_z, size_idx, openlem_nx, openlem_ny);
      //velocity_table.TableBase<dim,double>::reinit(size_idx);
      //// Indexes through x, y, and z.
      //if constexpr (dim == 2)
      //  {
      //    for (unsigned int x=0; x<velocity_table.size()[0]; ++x)
      //      for (unsigned int y=0; y<(velocity_table.size()[1]); ++y)
      //        // Convert back to m/s.
      //        velocity_table(x,y) = mesh_velocity_z[x][y];
      //  }
      //else
      //  {
      //    for (unsigned int x=0; x<velocity_table.size()[0]; ++x)
      //      for (unsigned int y=0; y<velocity_table.size()[1]; ++y)
      //        for (unsigned int z=0; z<velocity_table.size()[2]; ++z)
      //          // Convert back to m/s.
      //          velocity_table(x,y,z) = mesh_velocity_z[x][y];
      //  }
      //connector.updateXY();
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        connector.write_vtk(grid_extent[dim-1].second, this->get_timestep_number(),this->get_time()/year_in_seconds, dirname);
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 1)
        connector.write_vtk(grid_extent[dim-1].second, this->get_timestep_number(),this->get_time()/year_in_seconds, dirname, "mpi1");
      //std::cout << velocity_table.size()[0] << ":" << velocity_table.size()[1] << " : " << velocity_table.size()[2] << ", aplha = " << connector.alpha << ", dalpha = " << connector.dalpha << std::endl;
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years

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
      if (this->get_timestep_number() == 0)
        return;
      // Loop over all boundary indicators to set the velocity constraints
      /*for (const auto boundary_id : boundary_ids)
        VectorTools::interpolate_boundary_values (this->get_mapping(),
                                                  mesh_deformation_dof_handler,
                                                  boundary_id,
                                                  function,
                                                  mesh_velocity_constraints);
      */


      // As our grid_extent variable end points do not account for the change related to an origin
      // not at 0, we adjust this here into an interpolation extent.
      //std::array<std::pair<double,double>,dim> interpolation_extent;
      //for (unsigned int d=0; d<dim; ++d)
      //  {
      //    interpolation_extent[d].first = grid_extent[d].first;
      //    interpolation_extent[d].second = (grid_extent[d].second + grid_extent[d].first);
      //  }

      //Functions::InterpolatedUniformGridData<dim> *velocities;
      //Functions::InterpolatedUniformGridData<dim> velocities (interpolation_extent,
      //                                                        table_intervals,
      //                                                        velocity_table);

      VectorFunctionFromScalarFunctionObject<dim> vector_function_object(
        [&](const Point<dim> &p) -> double
      {
        //if (p[0] > -6345.87-100. && p[0] < -6345.87+15000. && p[1] > 28103. - 400. && p[1] < 28103. + 400. && p[2] > 25076.5 - 400. && p[2] < 25076.5+5000.)
        //  //if (p[0] >  3125.00-100. && p[0] <  3125.00+5000. && p[1] > 34375. - 400. && p[1] < 34375. + 400. && p[2] > 28028.4 - 400. && p[2] < 28028.4+5000.)
        //  {
        //    //std::cout << "p = " << p << ", vel = " << (connector.interpolation(p[0],p[1], aspect_mesh_elevation_h, aspect_mesh_z, aspect_dx, aspect_dy, grid_extent[0].first, grid_extent[1].first)-(p[2]-grid_extent[2].second)-(p[2]-grid_extent[2].second))/this->get_timestep() << ", interp = " << connector.interpolation(p[0],p[1], aspect_mesh_elevation_h, aspect_mesh_z, aspect_dx, aspect_dy, grid_extent[0].first, grid_extent[1].first)-(p[2]-grid_extent[2].second) << ", p2 = " << p[2] << ", ge = " << grid_extent[2].second << ", dt = " << this->get_timestep() << ", without timestep = " << (connector.interpolation(p[0],p[1], aspect_mesh_elevation_h, aspect_mesh_z, aspect_dx, aspect_dy, grid_extent[0].first, grid_extent[1].first)-(p[2]-grid_extent[2].second)-(p[2]-grid_extent[2].second)) << ", interp/timestep = " <<  (connector.interpolation(p[0],p[1], aspect_mesh_elevation_h, aspect_mesh_z, aspect_dx, aspect_dy, grid_extent[0].first, grid_extent[1].first))/this->get_timestep()<< std::endl;
        //  }
        //if ((connector.interpolation(p[0],p[1], aspect_mesh_elevation_h, aspect_mesh_z, aspect_dx, aspect_dy, grid_extent[0].first, grid_extent[1].first)+ grid_extent[dim-1].second )< 0.0)
        //if ((connector.interpolation(p[0],p[1], aspect_mesh_dh, aspect_mesh_z, aspect_dx, aspect_dy, grid_extent[0].first, grid_extent[1].first) )> 0.0)
        //  {
        //    std::cout << "rank = " <<  Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) << "p = " << p << ", vel = " << (connector.interpolation(p[0],p[1], aspect_mesh_dh, aspect_mesh_z, aspect_dx, aspect_dy, grid_extent[0].first, grid_extent[1].first))/this->get_timestep() << ", interp = " << connector.interpolation(p[0],p[1], aspect_mesh_dh, aspect_mesh_z, aspect_dx, aspect_dy, grid_extent[0].first, grid_extent[1].first) << ", p2 = " << p[2] << ", ge = " << grid_extent[2].second << ", dt = " << this->get_timestep() << ", without timestep = " << (connector.interpolation(p[0],p[1], aspect_mesh_dh, aspect_mesh_z, aspect_dx, aspect_dy, grid_extent[0].first, grid_extent[1].first)-(p[2]-grid_extent[2].second)-(p[2]-grid_extent[2].second)) << ", interp/timestep = " <<  (connector.interpolation(p[0],p[1], aspect_mesh_dh, aspect_mesh_z, aspect_dx, aspect_dy, grid_extent[0].first, grid_extent[1].first))/this->get_timestep()<< std::endl;
        //    //    AssertThrow(false, ExcMessage("wrong!!!!!!!"))  ;
        //  }
        return (connector.interpolation(p[0],p[1], aspect_mesh_dh, aspect_mesh_z, aspect_dx, aspect_dy, grid_extent[0].first, grid_extent[1].first))/this->get_timestep();//+ grid_extent[dim-1].second )/this->get_timestep();//-(p[2]-grid_extent[2].second))/this->get_timestep();// velocities.value(p);
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

      /*
          const types::boundary_id relevant_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id ("top");
          std::vector<std::vector<double>> local_aspect_values(dim+3, std::vector<double>());

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
                            const double indy = (vertex(1) - grid_extent[1].first)/openlem_dy;

                            if (std::abs(indy - std::round(indy)) >= node_tolerance)
                              continue;

                            //const double index = std::round((indy-1))*openlem_nx+std::round(indx);

                            local_aspect_values[0].push_back(vertex(dim-1) - grid_extent[dim-1].second);   //z component
                            if (indx == 5 && indy == 7)
                              std::cout << "vertex(dim-1) = " << vertex(dim-1) << ", ge[dim-1] = " << grid_extent[dim-1].second << std::endl;
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
      */
      return std::vector<std::vector<double>>();
    }
    template <int dim>
    void OpenLEM<dim>::execute_openlem(openlem::OceanGrid &grid,
                                       //std::vector<double> &elevation,
                                       //std::vector<double> &extra_vtk_field,
                                       //std::vector<double> &velocity_x,
                                       //std::vector<double> &velocity_y,
                                       //std::vector<double> &velocity_z,
                                       const double &openlem_timestep_in_years,
                                       const unsigned int &openlem_iterations,
                                       const std::string &dirname)
    {
      TimerOutput::Scope timer_section(this->get_computing_timer(), "Execute openlem");
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) != 0)
        return;

      // Because on the first timestep we will create an initial VTK file before running openlem
      // and a second after, we first set the visualization step to zero.
      unsigned int visualization_step = 0;
      const unsigned int current_timestep = this->get_timestep_number ();

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
        constexpr double sea_level = 0;

        for (unsigned int openlem_iteration = 0; openlem_iteration < openlem_iterations; ++openlem_iteration)//openlem_iterations; ++openlem_iteration)
          {
            //std::cout << "before = " << (grid.getNode(5,7))->h << std::endl;
            grid.computeWaterLevel();
            grid.clearOceans();
            grid.markOcean(deepest_point, sea_level);
            int nc = grid.computeFlowDirection();
            double ch = grid.erode(openlem_timestep_in_years);
            grid.findDeltas(openlem_timestep_in_years);
            grid.diffuse(openlem_timestep_in_years);
            printf("it: %i, Changes in flow direction: %i, Maximum elevation change: %e; ",openlem_iteration,nc,ch);
            unsigned int aspect_timestep_number = this->get_timestep_number();
            connector.updateCoordinateSystem (openlem_timestep_in_years);
            //connector.write_vtk(grid_extent[dim-1].second, aspect_timestep_number*10000+openlem_iteration,dirname);
            //std::cout << "after = " << (grid.getNode(5,7))->h << std::endl;
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
        std::cout << "deepest point = " << deepest_point.i << ":" << deepest_point.j << ", deepest point q = " <<  grid_new[deepest_point].q << ", deepest point h = " << grid_new[deepest_point].h  << ", b = " << (int)grid_new[deepest_point].b<< std::endl;
        //updateCoordinateSystem
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
    void OpenLEM<dim>::fill_openlem_arrays(openlem::OceanGrid &grid,
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
          openlem::Node *node = grid.getNode(index_x,index_y);
          node->h = local_aspect_values[0][i];
          node->l = 0;
          node->u = 2000;//1;//local_aspect_values[dim+2][i];
          //if (index_x == 5 && index_y == 7)
          //  std::cout << "Flag fill 1 before = " << node->h << ", " << local_aspect_values[0][i] << ", node_u = " << node->u << std::endl;
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
              openlem::Node *node = grid.getNode(index_x,index_y);
              node->h = local_aspect_values[0][i];
              node->l = 0;
              node->u = 1000;//local_aspect_values[dim+2][i];
              //if (index_x == 5 && index_y == 7)
              //  std::cout << "Flag fill 2 before = " << node->h << ", " << local_aspect_values[0][i] << ", node_u = " << node->u << std::endl;
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
        prm.enter_subsection ("OpenLEM");
        {
          prm.declare_entry("Minimize advection", "true",
                            Patterns::Bool(),
                            "Whether to rotate and translate the openLEM grid to minimize the advection");
          prm.declare_entry("Number of openlem timesteps per aspect timestep", "5",
                            Patterns::Integer(),
                            "Initial number of openlem time steps per ASPECT timestep, this value will double if"
                            "the openlem timestep is above the maximum openlem timestep.");
          prm.declare_entry("Kd", "1e-14",
                            Patterns::Double(0),
                            "Kd erosion rate for openlem. Units: $\\{yrs}$");
          prm.declare_entry("Kt", "1e-14",
                            Patterns::Double(0),
                            "Kt erosion rate for openlem. Units: $\\{yrs}$");
          prm.declare_entry("Ocean diffusivity", "1",
                            Patterns::Double(0),
                            "The diffusivisty of the ocean openlem. Units: $\\{yrs}$");
          prm.declare_entry("Maximum timestep length", "10e3",
                            Patterns::Double(0),
                            "Maximum timestep for openlem. Units: $\\{yrs}$");
          prm.declare_entry("Vertical exaggeration", "-1",
                            Patterns::Double(),
                            "Vertical exaggeration for openlem's VTK file. -1 outputs topography, basement, and sealevel.");
          prm.declare_entry("Additional openlem refinement", "0",
                            Patterns::Integer(),
                            "How many levels above ASPECT openlem should be refined.");
          prm.declare_entry ("Average out of plane surface topography in 2d", "true",
                             Patterns::Bool (),
                             "If this is set to false, then a 2D model will only consider the "
                             "center slice openlem gives. If set to true, then ASPECT will"
                             "average the mesh along Y excluding the ghost nodes.");
          prm.declare_entry("openlem seed", "1000",
                            Patterns::Integer(),
                            "Seed used for adding an initial noise to openlem topography based on the initial noise magnitude.");
          prm.declare_entry("Maximum surface refinement level", "1",
                            Patterns::Integer(),
                            "This should be set to the highest ASPECT refinement level expected at the surface.");
          prm.declare_entry("Surface refinement difference", "0",
                            Patterns::Integer(),
                            "The difference between the lowest and highest refinement level at the surface. E.g., if three resolution "
                            "levels are expected, this would be set to 2.");
          prm.declare_entry ("Use marine component", "false",
                             Patterns::Bool (),
                             "Flag to use the marine component of openlem.");
          prm.declare_entry("Y extent in 2d", "100000",
                            Patterns::Double(),
                            "openlem Y extent when using a 2D ASPECT model. Units: $\\{m}$");
          prm.declare_entry ("Use ghost nodes", "true",
                             Patterns::Bool (),
                             "Flag to use ghost nodes");
          prm.declare_entry ("Uplift and advect with openlem", "true",
                             Patterns::Bool (),
                             "Flag to use openlem advection and uplift.");
          prm.declare_entry("Node tolerance", "0.001",
                            Patterns::Double(),
                            "Node tolerance for how close an ASPECT node must be to the openlem node for the value to be transferred.");
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
            prm.declare_entry ("R ight", "1",
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
                               "Whether to set the ghost nodes at the openlem back and front boundary "
                               "to periodic even if 'Back' and 'Front' are set to fixed boundary.");
            prm.declare_entry ("Left right ghost nodes periodic", "false",
                               Patterns::Bool (),
                               "Whether to set the ghost nodes at the openlem left and right boundary "
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
                               "true will set all ghost nodes of fixed openlem boundaries to the height you specify in "
                               "'set Erosional base level'. \nThis can make "
                               "sense for a continental model where the model surrounding topography is assumed above sea level, "
                               "e.g. highlands. If the sea level would be used as an erosional base level in this case, all topography "
                               "erodes away with lots of 'sediment volume' lost through the sides of the model. This is mostly "
                               "important, when there are mountains in the middle of the model, while it is less important when there "
                               "is lower relief in the middle of the model. \n"
                               "In the openlem  visualization files, setting the extra base level may show up as a strong "
                               "slope at the fixed boundaries of the model. However, in the ASPECT visualization files it will not "
                               "show up, as the ghost nodes only exist in openlem.");
            prm.declare_entry("Erosional base level", "0",
                              Patterns::Double(),
                              "When 'Use a fixed erosional base level' is set to true, all ghost nodes of fixed "
                              "openlem boundaries where no mass flux is specified by the user (openlem boundary condition set to 1 "
                              "and 'Left/Right/Bottom/Top mass flux' set to 0) will be fixed to this elevation. The "
                              "reflecting boundaries (openlem boundary condition set to 0) will not be affected, nor are the "
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
    void OpenLEM<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection("OpenLEM");
        {
          openlem_ocean_diffusivity = prm.get_double("Ocean diffusivity");
          openlem_minimize_advection = prm.get_bool("Minimize advection");
          openlem_kd = prm.get_double("Kd");
          openlem_kt = prm.get_double("Kt");
          openlem_steps_per_aspect_step = prm.get_integer("Number of openlem timesteps per aspect timestep");
          maximum_openlem_timestep = prm.get_double("Maximum timestep length");
          //vexp = prm.get_double("Vertical exaggeration");
          additional_refinement_levels = prm.get_integer("Additional openlem refinement");
          //average_out_of_plane_surface_topography = prm.get_bool("Average out of plane surface topography in 2d");
          //openlem_seed = prm.get_integer("openlem seed");
          maximum_surface_refinement_level = prm.get_integer("Maximum surface refinement level");
          surface_refinement_difference = prm.get_integer("Surface refinement difference");
          //use_marine_component = prm.get_bool("Use marine component");
          openlem_y_extent_2d = prm.get_double("Y extent in 2d");
          use_ghost_nodes = prm.get_bool("Use ghost nodes");
          //openlem_advection_uplift = prm.get_bool("Uplift and advect with openlem");
          node_tolerance = prm.get_double("Node tolerance");
          //noise_elevation = prm.get_double("Initial noise magnitude");
          //sediment_rain_rates = Utilities::string_to_double
          //                      (Utilities::split_string_list(prm.get ("Sediment rain rates")));
          //sediment_rain_times = Utilities::string_to_double
          //                      (Utilities::split_string_list(prm.get ("Sediment rain time intervals")));

          //if (!this->convert_output_to_years())
          //  {
          //    maximum_openlem_timestep /= year_in_seconds;
          //    for (unsigned int j=0; j<sediment_rain_rates.size(); ++j)
          //      sediment_rain_rates[j] *= year_in_seconds;
          //  }

          //if (sediment_rain_rates.size() != sediment_rain_times.size()+1)
          //  AssertThrow(false, ExcMessage("Error: There must be one more sediment rain rate than time interval."));

          //for (unsigned int i=1; i<sediment_rain_times.size(); ++i)
          //  AssertThrow(sediment_rain_times[i] > sediment_rain_times[i-1], ExcMessage("Sediment rain time intervals must be an increasing array."));

          prm.enter_subsection("Boundary conditions");
          {
            // bottom = prm.get_integer("Front");
            // right = prm.get_integer("Right");
            // top = prm.get_integer("Back");
            // left = prm.get_integer("Left");
            // left_flux = prm.get_double("Left mass flux");
            // right_flux = prm.get_double("Right mass flux");
            // top_flux = prm.get_double("Back mass flux");
            // bottom_flux = prm.get_double("Front mass flux");

            // if (!this->convert_output_to_years())
            //   {
            //     left_flux *= year_in_seconds;
            //     right_flux *= year_in_seconds;
            //     top_flux *= year_in_seconds;
            //     bottom_flux *= year_in_seconds;
            //   }

            // // Put the boundary condition values into a four digit value to send to openlem.
            // openlem_boundary_conditions = bottom*1000+right*100+top*10+left;

            // if ((left_flux != 0 && top_flux != 0) || (left_flux != 0 && bottom_flux != 0) ||
            //     (right_flux != 0 && bottom_flux != 0) || (right_flux != 0 && top_flux != 0))
            //   AssertThrow(false,ExcMessage("Currently the plugin does not support mass flux through adjacent boundaries."));

            // topbottom_ghost_nodes_periodic = prm.get_bool("Back front ghost nodes periodic");
            // leftright_ghost_nodes_periodic = prm.get_bool("Left right ghost nodes periodic");
          }
          prm.leave_subsection();

          prm.enter_subsection("Erosional parameters");
          {
            //drainage_area_exponent_m = prm.get_double("Drainage area exponent");
            //slope_exponent_n = prm.get_double("Slope exponent");
            //sediment_river_incision_rate = prm.get_double("Sediment river incision rate");
            //bedrock_river_incision_rate = prm.get_double("Bedrock river incision rate");
            //sediment_transport_coefficient = prm.get_double("Sediment diffusivity");
            //bedrock_transport_coefficient = prm.get_double("Bedrock diffusivity");
            //bedrock_deposition_g = prm.get_double("Bedrock deposition coefficient");
            //sediment_deposition_g = prm.get_double("Sediment deposition coefficient");
            //slope_exponent_p = prm.get_double("Multi-direction slope exponent");
            //flat_elevation = prm.get_integer("Orographic elevation control");
            //wind_barrier_elevation = prm.get_integer("Orographic wind barrier height");
            //flat_erosional_factor = prm.get_double("Elevation factor");
            //wind_barrier_erosional_factor = prm.get_double("Wind barrier factor");
            //stack_controls = prm.get_bool("Stack orographic controls");
            //use_orographic_controls = prm.get_bool("Flag to use orographic controls");

            //if (!this->convert_output_to_years())
            //  {
            //    bedrock_river_incision_rate *= year_in_seconds;
            //    bedrock_transport_coefficient *= year_in_seconds;
            //    sediment_river_incision_rate *= year_in_seconds;
            //    bedrock_transport_coefficient *= year_in_seconds;
            //  }

            //// Wind direction
            //if (prm.get ("Wind direction") == "west")
            //  wind_direction = 0;
            //else if (prm.get ("Wind direction") == "east")
            //  wind_direction = 1;
            //else if (prm.get ("Wind direction") == "north")
            //  wind_direction = 2;
            //else if (prm.get ("Wind direction") == "south")
            //  wind_direction = 3;
            //else
            //  AssertThrow(false, ExcMessage("Not a valid wind direction."));

            //// set fixed ghost nodes to a base level for erosion that differs from sea level
            //use_fixed_erosional_base = prm.get_bool("Use a fixed erosional base level");
            //if (use_fixed_erosional_base)
            //  AssertThrow(use_fixed_erosional_base && use_ghost_nodes, ExcMessage(
            //                "If you want to use an erosional base level differing from sea level, "
            //                "you need to use ghost nodes."));
            //h_erosional_base = prm.get_double("Erosional base level");
          }
          prm.leave_subsection();

          prm.enter_subsection("Marine parameters");
          {
            //sea_level = prm.get_double("Sea level");
            //sand_surface_porosity = prm.get_double("Sand porosity");
            //silt_surface_porosity = prm.get_double("Silt porosity");
            //sand_efold_depth = prm.get_double("Sand e-folding depth");
            //silt_efold_depth = prm.get_double("Silt e-folding depth");
            //sand_silt_ratio = prm.get_double("Sand-silt ratio");
            //sand_silt_averaging_depth = prm.get_double("Depth averaging thickness");
            //sand_transport_coefficient = prm.get_double("Sand transport coefficient");
            //silt_transport_coefficient = prm.get_double("Silt transport coefficient");

            //if (!this->convert_output_to_years())
            //  {
            //    sand_transport_coefficient *= year_in_seconds;
            //    silt_transport_coefficient *= year_in_seconds;
            //  }
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
