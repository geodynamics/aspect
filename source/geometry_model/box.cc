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


#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/initial_topography_model/zero_topography.h>
#include <aspect/mesh_deformation/interface.h>
#include <aspect/simulator_signals.h>
#include <aspect/utilities.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/function_lib.h>


namespace aspect
{
  namespace GeometryModel
  {
    template <int dim>
    void
    Box<dim>::initialize ()
    {
      // Get pointer to initial topography model
      topo_model = const_cast<InitialTopographyModel::Interface<dim>*>(&this->get_initial_topography_model());

      // Check that initial topography is required.
      // If so, connect the initial topography function
      // to the right signals: It should be applied after
      // the final initial adaptive refinement and after a restart.
      if (Plugins::plugin_type_matches<InitialTopographyModel::ZeroTopography<dim>>(*topo_model) == false)
        {
          this->get_signals().pre_set_initial_state.connect(
            [&](typename parallel::distributed::Triangulation<dim> &tria)
          {
            this->topography(tria);
          }
          );
          this->get_signals().post_resume_load_user_data.connect(
            [&](typename parallel::distributed::Triangulation<dim> &tria)
          {
            this->topography(tria);
          }
          );
        }
    }


    template <int dim>
    void
    Box<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      const std::vector<unsigned int> rep_vec(repetitions.begin(), repetitions.end());
      GridGenerator::subdivided_hyper_rectangle (coarse_grid,
                                                 rep_vec,
                                                 box_origin,
                                                 box_origin+extents,
                                                 true);

      // Tell p4est about the periodicity of the mesh.
      std::vector<GridTools::PeriodicFacePair<typename parallel::distributed::Triangulation<dim>::cell_iterator>>
      periodicity_vector;
      for (int i=0; i<dim; ++i)
        if (periodic[i])
          GridTools::collect_periodic_faces
          ( coarse_grid, /*b_id1*/ 2*i, /*b_id2*/ 2*i+1,
            /*direction*/ i, periodicity_vector);

      if (periodicity_vector.size() > 0)
        coarse_grid.add_periodicity (periodicity_vector);
    }

    template <int dim>
    void
    Box<dim>::
    topography (typename parallel::distributed::Triangulation<dim> &grid) const
    {
      // Here we provide GridTools with the function to displace vertices
      // in the vertical direction by an amount specified by the initial topography model
      GridTools::transform(
        [&](const Point<dim> &p) -> Point<dim>
      {
        return this->add_topography(p);
      },
      grid);

      this->get_pcout() << "   Added initial topography to grid" << std::endl << std::endl;
    }


    template <int dim>
    Point<dim>
    Box<dim>::
    add_topography (const Point<dim> &x_y_z) const
    {
      // Get the surface x (,y) point
      Point<dim-1> surface_point;
      for (unsigned int d=0; d<dim-1; ++d)
        surface_point[d] = x_y_z[d];

      // Get the surface topography at this point
      const double topo = topo_model->value(surface_point);

      // Compute the displacement of the z coordinate
      const double ztopo = (x_y_z[dim-1] - box_origin[dim-1]) / extents[dim-1] * topo;

      // Compute the new point
      Point<dim> x_y_ztopo = x_y_z;
      x_y_ztopo[dim-1] += ztopo;

      return x_y_ztopo;
    }


    template <int dim>
    std::set<types::boundary_id>
    Box<dim>::
    get_used_boundary_indicators () const
    {
      // boundary indicators are zero through 2*dim-1
      std::set<types::boundary_id> s;
      for (unsigned int i=0; i<2*dim; ++i)
        s.insert (i);
      return s;
    }


    template <int dim>
    std::map<std::string,types::boundary_id>
    Box<dim>::
    get_symbolic_boundary_names_map () const
    {
      switch (dim)
        {
          case 2:
          {
            return
            {
              {"left",   0},
              {"right",  1},
              {"bottom", 2},
              {"top",    3}
            };
          }

          case 3:
          {
            return
            {
              {"left",   0},
              {"right",  1},
              {"front",  2},
              {"back",   3},
              {"bottom", 4},
              {"top",    5}

            };
          }
        }

      Assert (false, ExcNotImplemented());
      return {};
    }


    template <int dim>
    std::set<std::pair<std::pair<types::boundary_id, types::boundary_id>, unsigned int>>
    Box<dim>::
    get_periodic_boundary_pairs () const
    {
      std::set<std::pair<std::pair<types::boundary_id, types::boundary_id>, unsigned int>> periodic_boundaries;
      for ( unsigned int i=0; i<dim; ++i)
        if (periodic[i])
          periodic_boundaries.insert( std::make_pair( std::pair<types::boundary_id, types::boundary_id>(2*i, 2*i+1), i) );
      return periodic_boundaries;
    }



    template <int dim>
    void
    Box<dim>::adjust_positions_for_periodicity (Point<dim> &position,
                                                const ArrayView<Point<dim>> &connected_positions,
                                                const ArrayView<Tensor<1, dim>> &/*connected_velocities*/) const
    {
      for (unsigned int i = 0; i < dim; ++i)
        if (periodic[i])
          {
            if (position[i] < box_origin[i])
              {
                position[i] += extents[i];
                for (auto &connected_position: connected_positions)
                  connected_position[i] += extents[i];
              }
            else if (position[i] > box_origin[i] + extents[i])
              {
                position[i] -= extents[i];
                for (auto &connected_position: connected_positions)
                  connected_position[i] -= extents[i];
              }
          }
    }



    template <int dim>
    Point<dim>
    Box<dim>::get_extents () const
    {
      return extents;
    }

    template <int dim>
    const std::array<unsigned int, dim> &
    Box<dim>::get_repetitions () const
    {
      return repetitions;
    }

    template <int dim>
    Point<dim>
    Box<dim>::get_origin () const
    {
      return box_origin;
    }

    template <int dim>
    double
    Box<dim>::
    length_scale () const
    {
      return 0.01*extents[0];
    }


    template <int dim>
    double
    Box<dim>::depth(const Point<dim> &position) const
    {
      // Get the surface x (,y) point
      Point<dim-1> surface_point;
      for (unsigned int d=0; d<dim-1; ++d)
        surface_point[d] = position[d];

      // Get the surface topography at this point
      const double topo = topo_model->value(surface_point);

      const double d = extents[dim-1] + topo - (position(dim-1)-box_origin[dim-1]);
      return std::min (std::max (d, 0.), maximal_depth());
    }

    template <int dim>
    void
    Box<dim>::
    update ()
    {
        AssertThrow(Plugins::plugin_type_matches<const GeometryModel::Box<dim>>(this->get_geometry_model()),
                  ExcMessage("The geometry update function currently only works with a 2D box geometry model."));

        AssertThrow(dim==2,
                  ExcMessage("The geometry update function currently only works with a 2D box geometry model."));

        this->get_pcout() << "   Updating surface values... " << std::endl;
        TimerOutput::Scope timer_section(this->get_computing_timer(), "Geometry model surface update");

        // loop over all of the surface cells and save the elevation to a stored value.
        // This needs to be sent to 1 processor, sorted, and broadcast so that every processor knows the entire surface.
        // (Does bcast create memory of the entire variable on every processor or only a pointer to the one stored on 0?)
        std::vector<std::vector<double>> local_surface_height;
        Table<dim-1, double> data_table;
        const types::boundary_id relevant_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id ("top");
        const QTrapezoid<dim-1> face_corners;
        FEFaceValues<dim> fe_face_values(this->get_mapping(),
                                          this->get_fe(),
                                          face_corners,
                                          update_quadrature_points);

        // Loop over all corners at the surface and save their position.
        // TODO: Update this to work in 3D. Spherical?
        for (const auto &cell : this->get_dof_handler().active_cell_iterators())
          if (cell->is_locally_owned() && cell->at_boundary())
            for (const unsigned int face_no : cell->face_indices())
              if (cell->face(face_no)->at_boundary())
                {
                  if ( cell->face(face_no)->boundary_id() != relevant_boundary)
                    continue;

                  fe_face_values.reinit(cell, face_no);

                  for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
                    {
                      const Point<dim> vertex = fe_face_values.quadrature_point(corner);

                      // We can't push back a point so we convert it into a vector.
                      // This is needed later to keep the vertices together when sorting.
                      std::vector<double> vertex_row;
                      for(unsigned int i=0; i<dim; ++i)
                        vertex_row.push_back(vertex[i]);

                      local_surface_height.push_back(vertex_row);
                    }
                }

        // Combine all local_surfaces and broadcast back.
        std::vector<std::vector<double>> temp_surface;
        if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          // Save the surface stored on processor 0
          for (unsigned int i=0; i<local_surface_height.size(); i++)
                  temp_surface.push_back(local_surface_height[i]);

          // Save the surface stored on all other processors. will this need to be sorted
          // when using more than only a few processors?
          for (unsigned int p=1; p<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++p)
          {
            // First, find out the size of the array a process wants to send.
           MPI_Status status;
           int incoming_size = 0;
           MPI_Recv(&incoming_size, 1, MPI_INT, p, 42, this->get_mpi_communicator(), &status);

            std::vector<std::vector<double>> nonlocal_surface_height(incoming_size, std::vector<double>(dim));

            for (unsigned int i=0; i<nonlocal_surface_height.size(); ++i)
              MPI_Recv(&nonlocal_surface_height[i][0], dim, MPI_DOUBLE, p, 42, this->get_mpi_communicator(), &status);

            for (unsigned int i=0; i<nonlocal_surface_height.size(); ++i)
                    temp_surface.push_back(nonlocal_surface_height[i]);
          }

          // Sort the points. This may not be necessary in 2D.
          std::sort(temp_surface.begin(), temp_surface.end());

          // Define a comparison to remove duplicate surface points.
          bool (*compareRows)(const std::vector<double>&, const std::vector<double>&) = [](const std::vector<double>& row1, const std::vector<double>& row2) {
            return row1 == row2;
          };

          // Remove non-unique rows from the sorted 2D vector
          auto last = std::unique(temp_surface.begin(), temp_surface.end(), compareRows);
          temp_surface.erase(last, temp_surface.end());

          // First, bcast the size so that other processors can resize the relevant vector
          // to store the incoming data.
          int vector_size = temp_surface.size();

          TableIndices<dim-1> size_idx;
          for (unsigned int d=0; d<dim-1; ++d)
              size_idx[d] = vector_size;

          data_table.TableBase<dim-1,double>::reinit(size_idx);
          TableIndices<dim-1> idx;
			
          if(dim==2)
          {
              for (unsigned int x=0; x<(data_table.size()[0]); ++x)
              {
                idx[0] = x;
                data_table(idx) = temp_surface[x][1];
              }
          }

          Utilities::MPI::broadcast(this->get_mpi_communicator(), temp_surface, 0);
        }
        else
        {
          // Send local surface data.
          int vector_size = local_surface_height.size();
          MPI_Send(&vector_size, 1, MPI_INT, 0, 42, this->get_mpi_communicator());

          for (unsigned int i=0; i<local_surface_height.size(); i++)
              MPI_Ssend(&local_surface_height[i][0], dim, MPI_DOUBLE, 0, 42, this->get_mpi_communicator());

          // Do I need to resize the temp_surface variable when using this broadcast?
          temp_surface = Utilities::MPI::broadcast(this->get_mpi_communicator(), temp_surface, 0);
        }
		
        data_table.replicate_across_communicator (this->get_mpi_communicator(), 0);
		
		std::array<std::vector<double>, dim-1> coordinates;
        for(unsigned int i=0; i<temp_surface.size(); ++i)
        {
            coordinates[0].push_back(temp_surface[i][0]);
        }

		surface_function = new Functions::InterpolatedTensorProductGridData<dim-1> (coordinates, data_table);
    }

    template <int dim>
    double
    Box<dim>::depth_including_mesh_deformation(const Point<dim> &position) const
    {

      AssertThrow(Plugins::plugin_type_matches<const GeometryModel::Box<dim>>(this->get_geometry_model()),
                  ExcMessage("Depth with mesh deformation currently only works with a 2D box geometry model."));

      AssertThrow(dim==2,
                  ExcMessage("Depth with mesh deformation currently only works with a 2D box geometry model."));

      // Convert the point to dim-1, as we aren't interested in the vertical component
      // for the function.
      Point<dim-1> p;
      for (unsigned int d=0; d<dim-1; ++d)
        p[d] = position[d];

      double depth_from_surface = surface_function->value(p) - position[dim-1];
      return std::max (depth_from_surface, 0.);
    }


    template <int dim>
    double
    Box<dim>::height_above_reference_surface(const Point<dim> &position) const
    {
      return (position(dim-1)-box_origin[dim-1]) - extents[dim-1];
    }


    template <int dim>
    Point<dim>
    Box<dim>::representative_point(const double depth) const
    {
      Assert (depth >= 0,
              ExcMessage ("Given depth must be positive or zero."));
      Assert (depth <= maximal_depth(),
              ExcMessage ("Given depth must be less than or equal to the maximal depth of this geometry."));

      // choose a point on the center axis of the domain (without topography)
      Point<dim> p = extents/2+box_origin;

      // We need a dim-1 point to get the topo value.
      Point<dim-1> surface_point;
      for (unsigned int d=0; d<dim-1; ++d)
        surface_point[d] = p[d];

      const double topo = topo_model->value(surface_point);
      p[dim-1] = extents[dim-1]+box_origin[dim-1]-depth+topo;

      return p;
    }


    template <int dim>
    double
    Box<dim>::maximal_depth() const
    {
      return extents[dim-1] + topo_model->max_topography();
    }

    template <int dim>
    bool
    Box<dim>::has_curved_elements() const
    {
      return false;
    }

    template <int dim>
    bool
    Box<dim>::point_is_in_domain(const Point<dim> &point) const
    {
      // If mesh deformation is enabled, we have to loop over all the current
      // grid cells to see if the given point lies in the domain.
      // If mesh deformation is not enabled, or has not happened yet,
      // we can use the global extents of the model domain with or without
      // initial topography instead.
      // This function only checks if the given point lies in the domain
      // in its current shape at the current time. It can be called before
      // mesh deformation is applied in the first timestep (e.g., by the boundary
      // traction plugins), and therefore there is no guarantee
      // that the point will still lie in the domain after initial mesh deformation.
      if (this->get_parameters().mesh_deformation_enabled &&
          this->simulator_is_past_initialization())
        {
          return Utilities::point_is_in_triangulation<dim>(this->get_mapping(),
                                                           this->get_triangulation(),
                                                           point,
                                                           this->get_mpi_communicator());
        }
      // Without mesh deformation enabled, it is much cheaper to check whether the point lies in the domain.
      else
        {
          // The maximal extents of the unperturbed box domain.
          Point<dim> max_point = extents+box_origin;

          // If mesh deformation is not enabled, but initial topography
          // was/will be applied to the mesh, include this topography in the
          // extent of the domain.
          if (!Plugins::plugin_type_matches<const InitialTopographyModel::ZeroTopography<dim>>(this->get_initial_topography_model()))
            {
              // Get the surface x (,y) point
              Point<dim-1> surface_point;
              for (unsigned int d=0; d<dim-1; ++d)
                surface_point[d] = point[d];

              // Get the surface topography at this point
              const double topo = topo_model->value(surface_point);
              max_point[dim-1] += topo;
            }

          // Check whether point lies within the min/max coordinates of the domain including initial topography.
          for (unsigned int d = 0; d < dim; ++d)
            if (point[d] > max_point[d]+std::numeric_limits<double>::epsilon()*extents[d] ||
                point[d] < box_origin[d]-std::numeric_limits<double>::epsilon()*extents[d])
              return false;

          return true;
        }
    }

    template <int dim>
    std::array<double,dim>
    Box<dim>::cartesian_to_natural_coordinates(const Point<dim> &position_point) const
    {
      std::array<double,dim> position_array;
      for (unsigned int i = 0; i < dim; ++i)
        position_array[i] = position_point(i);

      return position_array;
    }


    template <int dim>
    aspect::Utilities::Coordinates::CoordinateSystem
    Box<dim>::natural_coordinate_system() const
    {
      return aspect::Utilities::Coordinates::CoordinateSystem::cartesian;
    }


    template <int dim>
    Point<dim>
    Box<dim>::natural_to_cartesian_coordinates(const std::array<double,dim> &position_tensor) const
    {
      Point<dim> position_point;
      for (unsigned int i = 0; i < dim; ++i)
        position_point[i] = position_tensor[i];

      return position_point;
    }


    template <int dim>
    void
    Box<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Box");
        {
          prm.declare_entry ("X extent", "1.",
                             Patterns::Double (0.),
                             "Extent of the box in x-direction. Units: \\si{\\meter}.");
          prm.declare_entry ("Y extent", "1.",
                             Patterns::Double (0.),
                             "Extent of the box in y-direction. Units: \\si{\\meter}.");
          prm.declare_entry ("Z extent", "1.",
                             Patterns::Double (0.),
                             "Extent of the box in z-direction. This value is ignored "
                             "if the simulation is in 2d. Units: \\si{\\meter}.");

          prm.declare_entry ("Box origin X coordinate", "0.",
                             Patterns::Double (),
                             "X coordinate of box origin. Units: \\si{\\meter}.");
          prm.declare_entry ("Box origin Y coordinate", "0.",
                             Patterns::Double (),
                             "Y coordinate of box origin. Units: \\si{\\meter}.");
          prm.declare_entry ("Box origin Z coordinate", "0.",
                             Patterns::Double (),
                             "Z coordinate of box origin. This value is ignored "
                             "if the simulation is in 2d. Units: \\si{\\meter}.");

          prm.declare_entry ("X repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in X direction.");
          prm.declare_entry ("Y repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in Y direction.");
          prm.declare_entry ("Z repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in Z direction.");

          prm.declare_entry ("X periodic", "false",
                             Patterns::Bool (),
                             "Whether the box should be periodic in X direction");
          prm.declare_entry ("Y periodic", "false",
                             Patterns::Bool (),
                             "Whether the box should be periodic in Y direction");
          prm.declare_entry ("Z periodic", "false",
                             Patterns::Bool (),
                             "Whether the box should be periodic in Z direction");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Box<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Box");
        {
          box_origin[0] = prm.get_double ("Box origin X coordinate");
          extents[0] = prm.get_double ("X extent");
          periodic[0] = prm.get_bool ("X periodic");
          repetitions[0] = prm.get_integer ("X repetitions");

          if (dim >= 2)
            {
              box_origin[1] = prm.get_double ("Box origin Y coordinate");
              extents[1] = prm.get_double ("Y extent");
              periodic[1] = prm.get_bool ("Y periodic");
              repetitions[1] = prm.get_integer ("Y repetitions");
            }

          if (dim >= 3)
            {
              // Use dim-1 instead of 2 to avoid compiler warning in 2d:
              box_origin[dim-1] = prm.get_double ("Box origin Z coordinate");
              extents[dim-1] = prm.get_double ("Z extent");
              periodic[dim-1] = prm.get_bool ("Z periodic");
              repetitions[dim-1] = prm.get_integer ("Z repetitions");
            }
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}


// explicit instantiations
namespace aspect
{
  namespace GeometryModel
  {
    ASPECT_REGISTER_GEOMETRY_MODEL(Box,
                                   "box",
                                   "A box geometry parallel to the coordinate directions. "
                                   "The extent of the box in each coordinate direction "
                                   "is set in the parameter file. The box geometry labels its "
                                   "2*dim sides as follows: in 2d, boundary indicators 0 through 3 "
                                   "denote the left, right, bottom and top boundaries; in 3d, boundary "
                                   "indicators 0 through 5 indicate left, right, front, back, bottom "
                                   "and top boundaries (see also the documentation of the deal.II class "
                                   "``ReferenceCell''). You can also use symbolic names ``left'', ``right'', "
                                   "etc., to refer to these boundaries in input files. "
                                   "It is also possible to add initial topography to the box model. Note however that "
                                   "this is done after the last initial adaptive refinement cycle. "
                                   "Also, initial topography is supposed to be small, as it is not taken into account "
                                   "when depth or a representative point is computed. ")


  }
}
