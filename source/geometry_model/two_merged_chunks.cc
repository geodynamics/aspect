/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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


#include <aspect/geometry_model/two_merged_chunks.h>
#include <aspect/geometry_model/initial_topography_model/zero_topography.h>
#include <aspect/geometry_model/initial_topography_model/ascii_data.h>

#include <aspect/simulator_signals.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>


namespace aspect
{
  namespace GeometryModel
  {
    template <int dim>
    void
    TwoMergedChunks<dim>::initialize ()
    {
      AssertThrow(Plugins::plugin_type_matches<const InitialTopographyModel::ZeroTopography<dim>>(this->get_initial_topography_model()) ||
                  Plugins::plugin_type_matches<const InitialTopographyModel::AsciiData<dim>>(this->get_initial_topography_model()),
                  ExcMessage("At the moment, only the Zero or AsciiData initial topography model can be used with the TwoMergedChunks geometry model."));

      manifold = std::make_unique<internal::ChunkGeometry<dim>>(this->get_initial_topography_model(),
                                                                 point1[1],
                                                                 point1[0],
                                                                 point2[0]-point1[0]);
    }



    template <int dim>
    void
    TwoMergedChunks<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      if (use_merged_grids)
        {
          // The two triangulations that will be merged into coarse_grid.
          Triangulation<dim> lower_coarse_grid;
          Triangulation<dim> upper_coarse_grid;

          const std::vector<unsigned int> lower_rep_vec(lower_repetitions.begin(), lower_repetitions.end());
          const std::vector<unsigned int> upper_rep_vec(upper_repetitions.begin(), upper_repetitions.end());

          // Create the lower box.
          GridGenerator::subdivided_hyper_rectangle (lower_coarse_grid,
                                                     lower_rep_vec,
                                                     point1,
                                                     point4,
                                                     false);

          // Create the upper box.
          GridGenerator::subdivided_hyper_rectangle (upper_coarse_grid,
                                                     upper_rep_vec,
                                                     point3,
                                                     point2,
                                                     false);

          // Merge the lower and upper mesh into one coarse_grid.
          // Now we have at least two cells.
          GridGenerator::merge_triangulations(lower_coarse_grid,
                                              upper_coarse_grid,
                                              coarse_grid);
        }
      else
        {
          const std::vector<unsigned int> lower_rep_vec(lower_repetitions.begin(), lower_repetitions.end());
          GridGenerator::subdivided_hyper_rectangle (coarse_grid,
                                                     lower_rep_vec,
                                                     point1,
                                                     point2,
                                                     false);

        }

      // Transform box into spherical chunk
      GridTools::transform (
        [&](const Point<dim> &p) -> Point<dim>
      {
        return manifold->push_forward(p);
      },
      coarse_grid);

      // Deal with a curved mesh by assigning a manifold. We arbitrarily
      // choose manifold_id 15 for this.
      coarse_grid.set_manifold (my_manifold_id, *manifold);
      for (const auto &cell : coarse_grid.active_cell_iterators())
        cell->set_all_manifold_ids (my_manifold_id);

      // Set the boundary indicators.
      set_boundary_indicators(coarse_grid);
      // Make sure the right boundary indicators are set after refinement
      // through the function set_boundary_indicators above.
      coarse_grid.signals.post_refinement.connect (
        [&]()
      {
        this->set_boundary_indicators(coarse_grid);
      });
    }



    template <int dim>
    void
    TwoMergedChunks<dim>::
    set_boundary_indicators (parallel::distributed::Triangulation<dim> &triangulation) const
    {
      // Iterate over all active cells and (re)set the boundary indicators.
      for (typename Triangulation<dim>::active_cell_iterator
           cell = triangulation.begin_active();
           cell != triangulation.end();
           ++cell)
        {
          // First set the default boundary indicators.
          for (const unsigned int f : cell->face_indices())
            if (cell->face(f)->at_boundary())
              cell->face(f)->set_boundary_id (f);

          if (cell->face(3)->at_boundary())
            // Set the upper part of the eastern boundary to indicator 2*dim+1
            // by comparing the radius of the face's center point with the
            // minimum radius of the upper radius.
            if (cell->face(3)->center(true).norm() > point3[0])
              cell->face(3)->set_boundary_id(2 * dim + 1);

          if (cell->face(2)->at_boundary())
            // Set the upper part of the western boundary to indicator 2*dim
            // by comparing the radius of the face's center point with the
            // minimum radius of the upper radius.
            if (cell->face(2)->center(true).norm() > point3[0])
              cell->face(2)->set_boundary_id (2*dim);

          if (dim==3)
            {
              // Set the upper part of the southern boundary to indicator 2*dim+2.
              if (cell->face(4)->at_boundary())
                if ((cell->vertex(cell->face(4)->n_vertices()-1).norm() + cell->vertex(0).norm()) / 2.0 > point3[0])
                  cell->face(4)->set_boundary_id (2*dim+2);

              // Set the upper part of the northern boundary to indicator 2*dim+3.
              if (cell->face(5)->at_boundary())
                if ((cell->vertex((cell->face(5)->n_vertices()-1)/2).norm()  + cell->vertex(0).norm()) / 2.0 > point3[0])
                  cell->face(5)->set_boundary_id (2*dim+3);
            }

        }
    }



    template <int dim>
    std::set<types::boundary_id>
    TwoMergedChunks<dim>::
    get_used_boundary_indicators () const
    {
      // boundary indicators are zero through 2*dim+2*(dim-1)-1
      std::set<types::boundary_id> s;
      for (unsigned int i=0; i<2*dim+2*(dim-1); ++i)
        s.insert (i);
      return s;
    }



    template <int dim>
    std::map<std::string,types::boundary_id>
    TwoMergedChunks<dim>::
    get_symbolic_boundary_names_map () const
    {
      switch (dim)
        {
          case 2:
          {
            return
            {
              {"bottom",     0},
              {"top",        1},
              {"lowerwest",  2},
              {"lowereast",  3},
              {"uppereast",  4},
              {"upperwest",  5}
            };
          }

          case 3:
          {
            return
            {
              {"bottom",     0},
              {"top",        1},
              {"lowerwest",  2},
              {"lowereast",  3},
              {"lowersouth", 4},
              {"lowernorth", 5},
              {"upperwest",  6},
              {"uppereast",  7},
              {"uppersouth", 8},
              {"uppernorth", 9}
            };
          }
        }

      Assert (false, ExcNotImplemented());
      return {};
    }



    template <int dim>
    double
    TwoMergedChunks<dim>::
    length_scale () const
    {
      // As described in the first ASPECT paper, a length scale of
      // 10km = 1e4m works well for the pressure scaling for earth
      // sized spherical shells. use a length scale that
      // yields this value for the R0,R1 corresponding to earth
      // but otherwise scales like (R1-R0)
      return 1e4 * maximal_depth() / (6336000.-3481000.);
    }



    template <int dim>
    double
    TwoMergedChunks<dim>::depth(const Point<dim> &position) const
    {
      // depth is defined wrt the reference surface point2[0]
      // negative depth is not allowed
      return std::max (0., std::min (point2[0]-position.norm(), maximal_depth()));
    }



    template <int dim>
    double
    TwoMergedChunks<dim>::height_above_reference_surface(const Point<dim> &position) const
    {
      return position.norm()-point2[0];
    }



    template <int dim>
    Point<dim>
    TwoMergedChunks<dim>::representative_point(const double depth) const
    {
      Assert (depth >= 0,
              ExcMessage ("Given depth must be positive or zero."));
      Assert (depth <= maximal_depth(),
              ExcMessage ("Given depth must be less than or equal to the maximal depth of this geometry."));

      // Choose a point at the mean longitude (and latitude)
      Point<dim> p = 0.5*(point2+point1);
      // at a depth beneath the top surface
      p[0] = point2[0]-depth;

      // Now convert to Cartesian coordinates
      return manifold->push_forward_sphere(p);
    }



    template <int dim>
    double
    TwoMergedChunks<dim>::west_longitude () const
    {
      return point1[1];
    }



    template <int dim>
    double
    TwoMergedChunks<dim>::east_longitude () const
    {
      return point2[1];
    }



    template <int dim>
    double
    TwoMergedChunks<dim>::longitude_range () const
    {
      return point2[1] - point1[1];
    }



    template <int dim>
    double
    TwoMergedChunks<dim>::south_latitude () const
    {
      if (dim == 3)
        return point1[2];
      else
        return 0;
    }



    template <int dim>
    double
    TwoMergedChunks<dim>::north_latitude () const
    {
      if (dim==3)
        return point2[2];
      else
        return 0;
    }



    template <int dim>
    double
    TwoMergedChunks<dim>::latitude_range () const
    {
      if (dim==3)
        return point2[2] - point1[2];
      else
        return 0;
    }



    template <int dim>
    double
    TwoMergedChunks<dim>::maximal_depth() const
    {
      return point2[0]-point1[0];
    }



    template <int dim>
    double
    TwoMergedChunks<dim>::inner_radius () const
    {
      return point1[0];
    }



    template <int dim>
    double
    TwoMergedChunks<dim>::outer_radius () const
    {
      return point2[0];
    }



    template <int dim>
    bool
    TwoMergedChunks<dim>::has_curved_elements() const
    {
      return true;
    }



    template <int dim>
    bool
    TwoMergedChunks<dim>::point_is_in_domain(const Point<dim> &point) const
    {
      AssertThrow(!this->get_parameters().mesh_deformation_enabled ||
                  // we are still before the first time step has started
                  this->get_timestep_number() == 0 ||
                  this->get_timestep_number() == numbers::invalid_unsigned_int,
                  ExcMessage("After displacement of the mesh, this function can no longer be used to determine whether a point lies in the domain or not."));

      AssertThrow(Plugins::plugin_type_matches<const InitialTopographyModel::ZeroTopography<dim>>(this->get_initial_topography_model()),
                  ExcMessage("After adding topography, this function can no longer be used to determine whether a point lies in the domain or not."));

      const Point<dim> spherical_point = manifold->pull_back(point);

      for (unsigned int d = 0; d < dim; ++d)
        if (spherical_point[d] > point2[d]+std::numeric_limits<double>::epsilon()*std::abs(point2[d]) ||
            spherical_point[d] < point1[d]-std::numeric_limits<double>::epsilon()*std::abs(point2[d]))
          return false;

      return true;
    }



    template <int dim>
    std::array<double,dim>
    TwoMergedChunks<dim>::cartesian_to_natural_coordinates(const Point<dim> &position_point) const
    {
      // the chunk manifold has a order of radius, longitude, latitude.
      // This is exactly what we need.
      // Ignore the topography to avoid a loop when calling the
      // AsciiDataBoundary for topography which uses this function....
      const Point<dim> transformed_point = manifold->pull_back_sphere(position_point);
      std::array<double,dim> position_array;
      for (unsigned int i = 0; i < dim; ++i)
        position_array[i] = transformed_point(i);

      return position_array;
    }



    template <int dim>
    aspect::Utilities::Coordinates::CoordinateSystem
    TwoMergedChunks<dim>::natural_coordinate_system() const
    {
      // TODO This will give problems somewhere down the line
      // if geometry models are asked for their coordinate system,
      // chunk returns spherical and then Utilities::Coordinates::cartesian_to_spherical
      // is used
      return aspect::Utilities::Coordinates::CoordinateSystem::spherical;
    }



    template <int dim>
    Point<dim>
    TwoMergedChunks<dim>::natural_to_cartesian_coordinates(const std::array<double,dim> &position_tensor) const
    {
      // Ignore the topography to avoid a loop when calling the
      // AsciiDataBoundary for topography which uses this function....
      Point<dim> position_point;
      for (unsigned int i = 0; i < dim; ++i)
        position_point[i] = position_tensor[i];
      const Point<dim> transformed_point = manifold->push_forward_sphere(position_point);

      return transformed_point;
    }



    template <int dim>
    void
    TwoMergedChunks<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Chunk with lithosphere boundary indicators");
        {
          prm.declare_entry ("Chunk inner radius", "0.",
                             Patterns::Double (0.),
                             "Radius at the bottom surface of the chunk. Units: \\si{\\meter}.");
          prm.declare_entry ("Chunk outer radius", "1.",
                             Patterns::Double (0.),
                             "Radius at the top surface of the chunk. Units: \\si{\\meter}.");
          prm.declare_entry ("Chunk middle boundary radius", "1",
                             Patterns::Double (0),
                             "Radius at the top surface of the lower chunk, "
                             "where it merges with the upper chunk. Units: \\si{\\meter}.");

          prm.declare_entry ("Chunk minimum longitude", "0.",
                             Patterns::Double (-180., 360.), // enables crossing of either hemisphere
                             "Minimum longitude of the chunk. Units: degrees.");
          prm.declare_entry ("Chunk maximum longitude", "1.",
                             Patterns::Double (-180., 360.), // enables crossing of either hemisphere
                             "Maximum longitude of the chunk. Units: degrees.");

          prm.declare_entry ("Chunk minimum latitude", "0.",
                             Patterns::Double (-90., 90.),
                             "Minimum latitude of the chunk. This value is ignored "
                             "if the simulation is in 2d. Units: degrees.");
          prm.declare_entry ("Chunk maximum latitude", "1.",
                             Patterns::Double (-90., 90.),
                             "Maximum latitude of the chunk. This value is ignored "
                             "if the simulation is in 2d. Units: degrees.");

          prm.declare_entry ("Outer chunk radius repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in radial direction for the upper chunk.");
          prm.declare_entry ("Inner chunk radius repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in radial direction for the lower chunk.");
          prm.declare_entry ("Longitude repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in longitude.");
          prm.declare_entry ("Latitude repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in latitude. This value is ignored "
                             "if the simulation is in 2d");

          prm.declare_entry ("Use merged grids", "true",
                             Patterns::Bool (),
                             "Whether to make the grid by gluing together two boxes, or just "
                             "use one chunk to make the grid. Using two grids glued together "
                             "is a safer option, since it forces the boundary conditions "
                             "to be always applied to the same depth, but using one grid allows "
                             "for a more flexible usage of the adaptive refinement. Note that if "
                             "there is no cell boundary exactly on the boundary between the lithosphere "
                             "and the mantle, the velocity boundary will not be exactly at that depth. "
                             "Therefore, using a merged "
                             "grid is generally recommended over using one grid. "
                             "When using one grid, the parameter for lower repetitions is used and the upper "
                             "repetitions are ignored.");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    TwoMergedChunks<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Chunk with lithosphere boundary indicators");
        {
          point1[0] = prm.get_double("Chunk inner radius");
          point2[0] = prm.get_double("Chunk outer radius");
          point3[0] = prm.get_double("Chunk middle boundary radius");
          point4[0] = point3[0];
          lower_repetitions[0] = prm.get_integer("Inner chunk radius repetitions");
          upper_repetitions[0] = prm.get_integer("Outer chunk radius repetitions");
          point1[1] = prm.get_double("Chunk minimum longitude") * constants::degree_to_radians;
          point2[1] = prm.get_double("Chunk maximum longitude") * constants::degree_to_radians;
          point3[1] = point1[1];
          point4[1] = point2[1];
          lower_repetitions[1] = prm.get_integer("Longitude repetitions");
          upper_repetitions[1] = lower_repetitions[1];

          AssertThrow(point1[0] < point2[0],
                      ExcMessage("Inner radius must be less than outer radius."));
          AssertThrow(point1[1] < point2[1],
                      ExcMessage("Minimum longitude must be less than maximum longitude."));
          AssertThrow(point2[1] - point1[1] < 2. * numbers::PI,
                      ExcMessage("Maximum - minimum longitude should be less than 360 degrees."));
          AssertThrow(point3[0] < point2[0],
                      ExcMessage("Middle boundary radius must be less than outer radius."));

          if (dim == 3)
            {
              point1[2] = prm.get_double ("Chunk minimum latitude") * constants::degree_to_radians;
              point2[2] = prm.get_double ("Chunk maximum latitude") * constants::degree_to_radians;
              point3[2] = point1[2];
              point4[2] = point2[2];
              lower_repetitions[2] = prm.get_integer ("Latitude repetitions");
              upper_repetitions[2] = lower_repetitions[2];

              AssertThrow (point1[2] < point2[2],
                           ExcMessage ("Minimum latitude must be less than maximum latitude."));
              AssertThrow (point1[2] > -0.5*numbers::PI,
                           ExcMessage ("Minimum latitude needs to be larger than -90 degrees."));
              AssertThrow (point2[2] < 0.5*numbers::PI,
                           ExcMessage ("Maximum latitude needs to be less than 90 degrees."));
            }
          use_merged_grids = prm.get_bool ("Use merged grids");
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
    ASPECT_REGISTER_GEOMETRY_MODEL(TwoMergedChunks,
                                   "chunk with lithosphere boundary indicators",
                                   "A geometry which can be described as a chunk of a spherical shell, "
                                   "bounded by lines of longitude, latitude and radius. "
                                   "The side boundaries have two boundary indicators, so the user "
                                   "can prescribe different boundary conditions on these boundaries. "
                                   "The minimum and maximum longitude, (latitude) and depth of the chunk "
                                   "are set in the parameter file. The chunk geometry labels its "
                                   "2*dim+2*(dim-1) sides as follows: ``lowerwest'' and ``lowereast'': "
                                   "minimum and maximum longitude of the lower part of the east and west "
                                   "side boundaries, ``upperwest'' and ``uppereast'': "
                                   "minimum and maximum longitude of the upper part of the east and west "
                                   "side boundaries, ``lowersouth'' and ``lowernorth'': "
                                   "minimum and maximum latitude of the lower part of the south and north "
                                   "side boundaries, ``uppersouth'' and ``uppernorth'': "
                                   "minimum and maximum latitude of the upper part of the south and north "
                                   "side boundaries, "
                                   "\n\n"
                                   "The dimensions of the model are specified by parameters "
                                   "of the following form: "
                                   "Chunk (minimum | maximum) (longitude | latitude): "
                                   "edges of geographical quadrangle (in degrees). "
                                   "Chunk (inner | outer | middle boundary) radius: Radii at bottom and top of chunk "
                                   "and the radius at which the lower boundary indicator along a side "
                                   "boundary transitions into the upper boundary indicator. "
                                   "(Longitude | Latitude) repetitions: "
                                   "number of cells in each coordinate direction."
                                   "(Inner | Outer) chunk radius repetitions: "
                                   "number of cells in the radial coordinate direction for the lower part "
                                   "of the domain (up to the Middle boundary radius) and for the upper part "
                                   "of the domain. "
                                   "\n\n"
                                   "When used in 2d, this geometry does not imply the use of "
                                   "a spherical coordinate system. Indeed, "
                                   "in 2d the geometry is simply a sector of an annulus in a Cartesian "
                                   "coordinate system and consequently would correspond to "
                                   "a sector of a cross section of the fluid filled space between two "
                                   "infinite cylinders where one has made the assumption that "
                                   "the velocity in direction of the cylinder axes is zero. "
                                   "This is consistent with the definition of what we consider "
                                   "the two-dimension case given in "
                                   "Section~\\ref{sec:methods:2d-models}. "
                                   "It is also possible to add initial topography to the chunk geometry, "
                                   "based on an ascii data file. ")
  }
}
