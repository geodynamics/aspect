/*
  Copyright (C) 2011 - 2017 by the authors of the ASPECT code.

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


#include <aspect/geometry_model/chunk.h>
#include <aspect/geometry_model/initial_topography_model/zero_topography.h>

#include <aspect/simulator_signals.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/manifold_lib.h>


namespace aspect
{
  namespace GeometryModel
  {
    template <int dim>
    Chunk<dim>::ChunkGeometry::ChunkGeometry()
      :
      point1_lon(0.0)
    {}

    template <int dim>
    DerivativeForm<1,dim,dim>
    Chunk<dim>::ChunkGeometry::
    push_forward_gradient(const Point<dim> &chart_point) const
    {
      const double R = chart_point[0]; // Radius
      const double phi = chart_point[1]; // Longitude

      Assert (R > 0.0, ExcMessage("Negative radius for given point."));

      DerivativeForm<1, dim, dim> DX;

      switch (dim)
        {
          case 2:
          {
            DX[0][0] =      std::cos(phi);
            DX[0][1] = -R * std::sin(phi);
            DX[1][0] =      std::sin(phi);
            DX[1][1] =  R * std::cos(phi);
            break;
          }
          case 3:
          {
            const double theta = chart_point[2]; // Latitude (not colatitude)

            DX[0][0] =      std::cos(theta) * std::cos(phi);
            DX[0][1] = -R * std::cos(theta) * std::sin(phi);
            DX[0][2] = -R * std::sin(theta) * std::cos(phi);
            DX[1][0] =      std::cos(theta) * std::sin(phi);
            DX[1][1] =  R * std::cos(theta) * std::cos(phi);
            DX[1][2] = -R * std::sin(theta) * std::sin(phi);
            DX[2][0] =      std::sin(theta);
            DX[2][1] = 0;
            DX[2][2] =  R * std::cos(theta);
            break;
          }
          default:
            Assert (false, ExcNotImplemented ());


        }

      return DX;
    }


    template <int dim>
    Point<dim>
    Chunk<dim>::ChunkGeometry::
    push_forward(const Point<dim> &input_vertex) const
    {
      Point<dim> output_vertex;
      switch (dim)
        {
          case 2:
          {
            output_vertex[0] = input_vertex[0]*std::cos(input_vertex[1]); // x=rcosphi
            output_vertex[1] = input_vertex[0]*std::sin(input_vertex[1]); // z=rsinphi
            break;
          }
          case 3:
          {
            output_vertex[0] = input_vertex[0]*std::cos(input_vertex[2])*std::cos(input_vertex[1]); // x=rsinthetacosphi
            output_vertex[1] = input_vertex[0]*std::cos(input_vertex[2])*std::sin(input_vertex[1]); // y=rsinthetasinphi
            output_vertex[2] = input_vertex[0]*std::sin(input_vertex[2]); // z=rcostheta
            break;
          }
          default:
            Assert (false, ExcNotImplemented ());
        }
      return output_vertex;
    }

    template <int dim>
    Point<dim>
    Chunk<dim>::ChunkGeometry::
    pull_back(const Point<dim> &v) const
    {
      Point<dim> output_vertex;
      switch (dim)
        {
          case 2:
          {
            output_vertex[1] = std::atan2(v[1], v[0]);
            output_vertex[0] = v.norm();
            // We must guarantee that all returned points have a longitude coordinate
            // value that is larger than (or equal to) the longitude of point1.
            // For example:
            // If the domain runs from longitude -10 to 200 degrees,
            // atan2 will also return a negative value (-180 to -160) for the points
            // with longitude 180 to 200. These values must be corrected
            // so that they are larger than the minimum longitude value of -10,
            // by adding 360 degrees.
            // A 100*epsilon ensures we catch all cases.
            if (output_vertex[1] < 0.0)
              if (output_vertex[1] < point1_lon - 100 * std::abs(point1_lon)*std::numeric_limits<double>::epsilon())
                output_vertex[1] += 2.0 * numbers::PI;
            break;
          }
          case 3:
          {
            const double radius=v.norm();
            output_vertex[0] = radius;
            output_vertex[1] = std::atan2(v[1], v[0]);
            // See 2D case
            if (output_vertex[1] < 0.0)
              if (output_vertex[1] < point1_lon - 100 * std::abs(point1_lon)*std::numeric_limits<double>::epsilon())
                output_vertex[1] += 2.0 * numbers::PI;
            output_vertex[2] = std::asin(v[2]/radius);
            break;
          }
          default:
            Assert (false, ExcNotImplemented ());
        }
      return output_vertex;
    }

    template <int dim>
    void
    Chunk<dim>::ChunkGeometry::
    set_min_longitude(const double p1_lon)
    {
      point1_lon = p1_lon;
    }

#if !DEAL_II_VERSION_GTE(9,0,0)
    template <int dim>
    void
    Chunk<dim>::initialize ()
    {
      // Call function to connect the set/clear manifold id functions
      // to the right signal
      connect_to_signal(this->get_signals());

    }
#endif

    template <int dim>
    void
    Chunk<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      std::vector<unsigned int> rep_vec(repetitions, repetitions+dim);
      GridGenerator::subdivided_hyper_rectangle (coarse_grid,
                                                 rep_vec,
                                                 point1,
                                                 point2,
                                                 true);

#if !DEAL_II_VERSION_GTE(9,0,0)
      // At this point, all boundary faces have their correct boundary
      // indicators, but the edges do not. We want the edges of curved
      // faces to be curved as well, so we set the edge boundary indicators
      // to the same boundary indicators as their faces.
      for (typename Triangulation<dim>::active_cell_iterator
           cell = coarse_grid.begin_active();
           cell != coarse_grid.end(); ++cell)
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->face(f)->at_boundary())
            // Set edges on the radial faces; both adjacent faces
            // should agree on where new points along the boundary lie
            // for these edges, so the order of the boundaries does not matter
            if ((cell->face(f)->boundary_id() == 2)
                ||
                (cell->face(f)->boundary_id() == 3))
              cell->face(f)->set_all_boundary_ids(cell->face(f)->boundary_id());

      for (typename Triangulation<dim>::active_cell_iterator
           cell = coarse_grid.begin_active();
           cell != coarse_grid.end(); ++cell)
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->face(f)->at_boundary())
            // Set edges on the radial faces; both adjacent faces
            // should agree on where new points along the boundary lie
            // for these edges, so the order of the boundaries does not matter
            if ((cell->face(f)->boundary_id() == 4)
                ||
                (cell->face(f)->boundary_id() == 5))
              cell->face(f)->set_all_boundary_ids(cell->face(f)->boundary_id());

      for (typename Triangulation<dim>::active_cell_iterator
           cell = coarse_grid.begin_active();
           cell != coarse_grid.end(); ++cell)
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->face(f)->at_boundary())
            // (Re-)Set edges on the spherical shells to ensure that
            // they are all curved as expected
            if ((cell->face(f)->boundary_id() == 0)
                ||
                (cell->face(f)->boundary_id() == 1))
              cell->face(f)->set_all_boundary_ids(cell->face(f)->boundary_id());
#endif

      // Transform box into spherical chunk
      GridTools::transform (std_cxx11::bind(&ChunkGeometry::push_forward,
                                            std_cxx11::cref(manifold),
                                            std_cxx11::_1),
                            coarse_grid);

      // Deal with a curved mesh
      // Attach the real manifold to slot 15.
      coarse_grid.set_manifold (15, manifold);
      for (typename Triangulation<dim>::active_cell_iterator cell =
             coarse_grid.begin_active(); cell != coarse_grid.end(); ++cell)
        cell->set_all_manifold_ids (15);

#if !DEAL_II_VERSION_GTE(9,0,0)
      // On the boundary faces, set boundary objects.
      // The east and west boundaries are straight,
      // the inner and outer boundary are part of
      // a spherical shell. In 3D, the north and south boundaries
      // are part of a cone with the tip at the origin.
      static const StraightBoundary<dim> boundary_straight;

      // Attach boundary objects to the straight east and west boundaries
      coarse_grid.set_boundary(2, boundary_straight);
      coarse_grid.set_boundary(3, boundary_straight);

      if (dim == 3)
        {
          // Define the center point of the greater radius end of the
          // north and south boundary cones.
          // These lie along the z-axis.
          Point<dim> center;
          Point<dim> north, south;
          const double outer_radius = point2[0];
          north[dim-1] = outer_radius * std::sin(point2[2]);
          south[dim-1] = outer_radius * std::sin(point1[2]);
          // Define the radius of the cones
          const double north_radius = std::sqrt(outer_radius*outer_radius-north[dim-1]*north[dim-1]);
          const double south_radius = std::sqrt(outer_radius*outer_radius-south[dim-1]*south[dim-1]);
          static const ConeBoundary<dim> boundary_cone_north(0.0,north_radius,center,north);
          static const ConeBoundary<dim> boundary_cone_south(0.0,south_radius,center,south);

          // Attach boundary objects to the conical north and south boundaries
          // If one of the boundaries lies at the equator,
          // just use the straight boundary.
          if (point2[2] != 0.0)
            coarse_grid.set_boundary (5, boundary_cone_north);
          else
            coarse_grid.set_boundary (5, boundary_straight);

          if (point1[2] != 0.0)
            coarse_grid.set_boundary (4, boundary_cone_south);
          else
            coarse_grid.set_boundary (4, boundary_straight);
        }

      // Attach shell boundary objects to the curved inner and outer boundaries
      static const HyperShellBoundary<dim> boundary_shell;
      coarse_grid.set_boundary (0, boundary_shell);
      coarse_grid.set_boundary (1, boundary_shell);
#endif
    }

#if !DEAL_II_VERSION_GTE(9,0,0)
    template <int dim>
    void
    Chunk<dim>::set_manifold_ids (typename parallel::distributed::Triangulation<dim> &triangulation)
    {
      // Set all cells, faces and edges to manifold_id 15
      for (typename Triangulation<dim>::active_cell_iterator cell =
             triangulation.begin_active(); cell != triangulation.end(); ++cell)
        cell->set_all_manifold_ids (15);
    }

    template <int dim>
    void
    Chunk<dim>::clear_manifold_ids (typename parallel::distributed::Triangulation<dim> &triangulation)
    {
      // Clear the manifold_id from the faces and edges at the boundary
      // so that the boundary objects can be used
      for (typename Triangulation<dim>::active_cell_iterator cell =
             triangulation.begin_active(); cell != triangulation.end(); ++cell)
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->face(f)->at_boundary())
            cell->face(f)->set_all_manifold_ids (numbers::invalid_manifold_id);
    }

    template <int dim>
    void
    Chunk<dim>::connect_to_signal (SimulatorSignals<dim> &signals)
    {
      // Connect the topography function to the signal
      signals.pre_compute_no_normal_flux_constraints.connect (std_cxx11::bind (&Chunk<dim>::clear_manifold_ids,
                                                                               std_cxx11::ref(*this),
                                                                               std_cxx11::_1));
      signals.post_compute_no_normal_flux_constraints.connect (std_cxx11::bind (&Chunk<dim>::set_manifold_ids,
                                                                                std_cxx11::ref(*this),
                                                                                std_cxx11::_1));
    }
#endif

    template <int dim>
    std::set<types::boundary_id>
    Chunk<dim>::
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
    Chunk<dim>::
    get_symbolic_boundary_names_map () const
    {
      switch (dim)
        {
          case 2:
          {
            static const std::pair<std::string,types::boundary_id> mapping[]
              = { std::pair<std::string,types::boundary_id>("bottom", 0),
                  std::pair<std::string,types::boundary_id>("top",    1),
                  std::pair<std::string,types::boundary_id>("west",   2),
                  std::pair<std::string,types::boundary_id>("east",   3)
                };

            return std::map<std::string,types::boundary_id> (&mapping[0],
                                                             &mapping[sizeof(mapping)/sizeof(mapping[0])]);
          }

          case 3:
          {
            static const std::pair<std::string,types::boundary_id> mapping[]
              = { std::pair<std::string,types::boundary_id>("bottom", 0),
                  std::pair<std::string,types::boundary_id>("top",    1),
                  std::pair<std::string,types::boundary_id>("west",   2),
                  std::pair<std::string,types::boundary_id>("east",   3),
                  std::pair<std::string,types::boundary_id>("south",  4),
                  std::pair<std::string,types::boundary_id>("north",  5)
                };

            return std::map<std::string,types::boundary_id> (&mapping[0],
                                                             &mapping[sizeof(mapping)/sizeof(mapping[0])]);
          }
        }

      Assert (false, ExcNotImplemented());
      return std::map<std::string,types::boundary_id>();
    }


    template <int dim>
    double
    Chunk<dim>::
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
    Chunk<dim>::depth(const Point<dim> &position) const
    {
      return std::min (std::max (point2[0]-position.norm(), 0.), maximal_depth());
    }


    template <int dim>
    Point<dim>
    Chunk<dim>::representative_point(const double depth) const
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
      return manifold.push_forward(p);
    }

    template <int dim>
    double
    Chunk<dim>::west_longitude () const
    {
      return point1[1];
    }


    template <int dim>
    double
    Chunk<dim>::east_longitude () const
    {
      return point2[1];
    }

    template <int dim>
    double
    Chunk<dim>::longitude_range () const
    {
      return point2[1] - point1[1];
    }

    template <int dim>
    double
    Chunk<dim>::south_latitude () const
    {
      if (dim == 3)
        return point1[2];
      else
        return 0;
    }


    template <int dim>
    double
    Chunk<dim>::north_latitude () const
    {
      if (dim==3)
        return point2[2];
      else
        return 0;
    }


    template <int dim>
    double
    Chunk<dim>::latitude_range () const
    {
      if (dim==3)
        return point2[2] - point1[2];
      else
        return 0;
    }


    template <int dim>
    double
    Chunk<dim>::maximal_depth() const
    {
      return point2[0]-point1[0];
    }

    template <int dim>
    double
    Chunk<dim>::inner_radius () const
    {
      return point1[0];
    }

    template <int dim>
    double
    Chunk<dim>::outer_radius () const
    {
      return point2[0];
    }

    template <int dim>
    bool
    Chunk<dim>::has_curved_elements() const
    {
      return true;
    }

    template <int dim>
    bool
    Chunk<dim>::point_is_in_domain(const Point<dim> &point) const
    {
      AssertThrow(this->get_free_surface_boundary_indicators().size() == 0 ||
                  this->get_timestep_number() == 0,
                  ExcMessage("After displacement of the free surface, this function can no longer be used to determine whether a point lies in the domain or not."));

      AssertThrow(dynamic_cast<const InitialTopographyModel::ZeroTopography<dim>*>(&this->get_initial_topography_model()) != 0,
                  ExcMessage("After adding topography, this function can no longer be used to determine whether a point lies in the domain or not."));

      const Point<dim> spherical_point = manifold.pull_back(point);

      for (unsigned int d = 0; d < dim; d++)
        if (spherical_point[d] > point2[d]+std::numeric_limits<double>::epsilon()*std::abs(point2[d]) ||
            spherical_point[d] < point1[d]-std::numeric_limits<double>::epsilon()*std::abs(point2[d]))
          return false;

      return true;
    }


    template <int dim>
    void
    Chunk<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Chunk");
        {
          prm.declare_entry ("Chunk inner radius", "0",
                             Patterns::Double (0),
                             "Radius at the bottom surface of the chunk. Units: m.");
          prm.declare_entry ("Chunk outer radius", "1",
                             Patterns::Double (0),
                             "Radius at the top surface of the chunk. Units: m.");

          prm.declare_entry ("Chunk minimum longitude", "0",
                             Patterns::Double (-180, 360), // enables crossing of either hemisphere
                             "Minimum longitude of the chunk. Units: degrees.");
          prm.declare_entry ("Chunk maximum longitude", "1",
                             Patterns::Double (-180, 360), // enables crossing of either hemisphere
                             "Maximum longitude of the chunk. Units: degrees.");

          prm.declare_entry ("Chunk minimum latitude", "0",
                             Patterns::Double (-90, 90),
                             "Minimum latitude of the chunk. This value is ignored "
                             "if the simulation is in 2d. Units: degrees.");
          prm.declare_entry ("Chunk maximum latitude", "1",
                             Patterns::Double (-90, 90),
                             "Maximum latitude of the chunk. This value is ignored "
                             "if the simulation is in 2d. Units: degrees.");

          prm.declare_entry ("Radius repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in radius.");
          prm.declare_entry ("Longitude repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in longitude.");
          prm.declare_entry ("Latitude repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in latitude. This value is ignored "
                             "if the simulation is in 2d");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Chunk<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Chunk");
        {

          const double degtorad = dealii::numbers::PI/180;

          Assert (dim >= 2, ExcInternalError());
          Assert (dim <= 3, ExcInternalError());

          if (dim >= 2)
            {
              point1[0] = prm.get_double ("Chunk inner radius");
              point2[0] = prm.get_double ("Chunk outer radius");
              repetitions[0] = prm.get_integer ("Radius repetitions");
              point1[1] = prm.get_double ("Chunk minimum longitude") * degtorad;
              point2[1] = prm.get_double ("Chunk maximum longitude") * degtorad;
              repetitions[1] = prm.get_integer ("Longitude repetitions");

              AssertThrow (point1[0] < point2[0],
                           ExcMessage ("Inner radius must be less than outer radius."));
              AssertThrow (point1[1] < point2[1],
                           ExcMessage ("Minimum longitude must be less than maximum longitude."));
              AssertThrow (point2[1] - point1[1] < 2.*numbers::PI,
                           ExcMessage ("Maximum - minimum longitude should be less than 360 degrees."));
            }

          // Inform the manifold about the minimum longitude
          manifold.set_min_longitude(point1[1]);

          if (dim == 3)
            {
              point1[2] = prm.get_double ("Chunk minimum latitude") * degtorad;
              point2[2] = prm.get_double ("Chunk maximum latitude") * degtorad;
              repetitions[2] = prm.get_integer ("Latitude repetitions");

              AssertThrow (point1[2] < point2[2],
                           ExcMessage ("Minimum latitude must be less than maximum latitude."));
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
    ASPECT_REGISTER_GEOMETRY_MODEL(Chunk,
                                   "chunk",
                                   "A geometry which can be described as a chunk of a spherical shell, "
                                   "bounded by lines of longitude, latitude and radius. "
                                   "The minimum and maximum longitude, (latitude) and depth of the chunk "
                                   "is set in the parameter file. The chunk geometry labels its "
                                   "2*dim sides as follows: ``west'' and ``east'': minimum and maximum "
                                   "longitude, ``south'' and ``north'': minimum and maximum latitude, "
                                   "``inner'' and ``outer'': minimum and maximum radii. "
                                   "Names in the parameter files are as follows: "
                                   "Chunk (minimum || maximum) (longitude || latitude): "
                                   "edges of geographical quadrangle (in degrees)"
                                   "Chunk (inner || outer) radius: Radii at bottom and top of chunk"
                                   "(Longitude || Latitude || Radius) repetitions: "
                                   "number of cells in each coordinate direction.")
  }
}

