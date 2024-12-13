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


#include <aspect/geometry_model/ellipsoidal_chunk.h>
#include <aspect/utilities.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <boost/lexical_cast.hpp>
#include <aspect/compat.h>


/**
 * This geometry model implements an (3d) ellipsoidal_chunk geometry which can be non-coordinate parallel.
 * @author This plugin is a joint effort of Menno Fraters, D Sarah Stamps and Wolfgang Bangerth
 */

namespace aspect
{
  namespace GeometryModel
  {
    namespace
    {
      template <int dim>
      void
      set_manifold_ids(Triangulation<dim> &triangulation)
      {
        for (const auto &cell : triangulation.active_cell_iterators())
          cell->set_all_manifold_ids (15);
      }

      template <int dim>
      void
      clear_manifold_ids(Triangulation<dim> &triangulation)
      {
        for (const auto &cell : triangulation.active_cell_iterators())
          cell->set_all_manifold_ids (numbers::flat_manifold_id);
      }
    }


    namespace internal
    {
      /*
       * the EllipsoidalChunkGeometry class
       */

      // Constructor
      template <int dim>
      EllipsoidalChunkGeometry<dim>::EllipsoidalChunkGeometry(const InitialTopographyModel::Interface<dim> &topo,
                                                              const double para_semi_major_axis_a,
                                                              const double para_eccentricity,
                                                              const double para_semi_minor_axis_b,
                                                              const double para_bottom_depth,
                                                              const std::vector<Point<2>> &para_corners)
        :
        topography (&topo),
        semi_major_axis_a (para_semi_major_axis_a),
        eccentricity (para_eccentricity),
        semi_minor_axis_b (para_semi_minor_axis_b),
        bottom_depth (para_bottom_depth),
        corners (para_corners)
      {
        AssertThrow (dim == 3, ExcMessage("This manifold can currently only be used in 3d."));
      }



      template <int dim>
      Point<3>
      EllipsoidalChunkGeometry<dim>::push_forward_ellipsoid(const Point<3> &phi_theta_d, const double semi_major_axis_a, const double eccentricity) const
      {
        // The following converts phi, theta and negative depth to x, y, z
        // Depth is measured perpendicular to the ellipsoid surface
        // (i.e. along a vector which does not generally pass through the origin)
        // Expressions can be found in Ellipsoidal and Cartesian Coordinates Conversion
        // Subirana, Zornoza and Hernandez-Pajares, 2011:
        // https://gssc.esa.int/navipedia/index.php/Ellipsoidal_and_Cartesian_Coordinates_Conversion

        const double phi   = phi_theta_d[0]; // Longitude in radians
        const double theta = phi_theta_d[1]; // Latitude in radians
        const double d     = phi_theta_d[2]; // The negative depth (a depth of 10 meters is -10)

        const double R_bar = semi_major_axis_a / std::sqrt(1 - (eccentricity * eccentricity *
                                                                std::sin(theta) * std::sin(theta))); // radius of curvature of the prime vertical

        return Point<3> ((R_bar + d) * std::cos(phi) * std::cos(theta),
                         (R_bar + d) * std::sin(phi) * std::cos(theta),
                         ((1 - eccentricity * eccentricity) * R_bar + d) * std::sin(theta));
      }

      template <int dim>
      Point<3>
      EllipsoidalChunkGeometry<dim>::pull_back_ellipsoid(const Point<3> &x, const double semi_major_axis_a, const double eccentricity) const
      {
        // The following converts x, y, z to phi, theta and negative depth
        // Depth is measured perpendicular to the ellipsoid surface
        // (i.e. along a vector which does not generally pass through the origin)
        // Expressions can be found in Ellipsoidal and Cartesian Coordinates Conversion
        // Subirana, Zornoza and Hernandez-Pajares, 2011:
        // https://gssc.esa.int/navipedia/index.php/Ellipsoidal_and_Cartesian_Coordinates_Conversion

        const double R      = semi_major_axis_a; // semi-major axis
        const double b      = R * std::sqrt(1 - eccentricity * eccentricity); // semi-minor axis
        const double p      = std::sqrt(x(0) * x(0) + x(1) * x(1)); // distance from origin projected onto x-y plane
        const double th     = std::atan2(R * x(2), b * p); // starting guess for theta
        const double phi    = std::atan2(x(1), x(0)); // azimuth (geodetic longitude)
        const double theta  = std::atan2(x(2) + (R * R - b * b) / b * Utilities::fixed_power<3>(std::sin(th)),
                                         (p - (eccentricity * eccentricity * R  * Utilities::fixed_power<3>(std::cos(th))))); // first iterate for theta
        const double R_bar  = R / (std::sqrt(1 - eccentricity * eccentricity * std::sin(theta) * std::sin(theta))); // first iterate for R_bar

        Point<3> phi_theta_d;
        phi_theta_d[0] = phi;

        phi_theta_d[1] = theta;
        phi_theta_d[2] = p / std::cos(theta) - R_bar; // first iterate for d
        return phi_theta_d;
      }

      template <int dim>
      Point<3>
      EllipsoidalChunkGeometry<dim>::push_forward_topography(const Point<3> &phi_theta_d_hat) const
      {
        const double d_hat = phi_theta_d_hat[2]; // long, lat, depth
        Point<dim-1> phi_theta;
        if (dim == 3)
          phi_theta = Point<dim-1>(phi_theta_d_hat[0] * constants::radians_to_degree,phi_theta_d_hat[1] * constants::radians_to_degree);
        const double h = topography != nullptr ? topography->value(phi_theta) : 0;
        const double d = d_hat + (d_hat + bottom_depth)/bottom_depth*h;
        const Point<3> phi_theta_d (phi_theta_d_hat[0],
                                    phi_theta_d_hat[1],
                                    d);
        return phi_theta_d;
      }

      template <int dim>
      Point<3>
      EllipsoidalChunkGeometry<dim>::pull_back_topography(const Point<3> &phi_theta_d) const
      {
        const double d = phi_theta_d[2];
        Point<dim-1> phi_theta;
        if (dim == 3)
          phi_theta = Point<dim-1>(phi_theta_d[0] * constants::radians_to_degree,
                                   phi_theta_d[1] * constants::radians_to_degree);
        const double h = topography != nullptr ? topography->value(phi_theta) : 0;
        const double d_hat = bottom_depth * (d-h)/(bottom_depth+h);
        const Point<3> phi_theta_d_hat (phi_theta_d[0],
                                        phi_theta_d[1],
                                        d_hat);
        return phi_theta_d_hat;
      }

      /**
       * TODO: These functions (pull back and push forward) should be changed that they always
       * take and return 3d points, because 2d points make no sense for an ellipsoid, even with
       * a 2d triangulation. To do this correctly we need to add the spacedim to the triangulation
       * in ASPECT. What is now presented is just a temporary fix to get access to the pull back
       * function from outside. The push forward function can't be fixed in this way, because
       * it is used by a bind statement.
       */
      template <int dim>
      Point<3>
      EllipsoidalChunkGeometry<dim>::pull_back(const Point<3> &space_point) const
      {
        return pull_back_topography(pull_back_ellipsoid (space_point, semi_major_axis_a, eccentricity));

      }

      template <int dim>
      Point<2>
      EllipsoidalChunkGeometry<dim>::pull_back(const Point<2> &space_point) const
      {
        return space_point;

      }

      template <int dim>
      Point<3>
      EllipsoidalChunkGeometry<dim>::push_forward(const Point<3> &chart_point) const
      {
        return push_forward_ellipsoid (push_forward_topography(chart_point), semi_major_axis_a, eccentricity);
      }

      template <int dim>
      std::unique_ptr<Manifold<dim,3>>
      EllipsoidalChunkGeometry<dim>::clone() const
      {
        return std::make_unique<EllipsoidalChunkGeometry>(*this);
      }
    }

    template <int dim>
    void
    EllipsoidalChunk<dim>::initialize()
    {
      manifold = std::make_unique<internal::EllipsoidalChunkGeometry<dim>>(this->get_initial_topography_model(),
                                                                            semi_major_axis_a,
                                                                            eccentricity,
                                                                            semi_minor_axis_b,
                                                                            bottom_depth,
                                                                            corners);
    }



    template <>
    void
    EllipsoidalChunk<3>::create_coarse_mesh(parallel::distributed::Triangulation<3> &coarse_grid) const
    {
      const int dim = 3;

      // Generate parallelepiped grid with one point (point 0) at (0,0,0) and the
      // other corners (respectively corner 1,2 and 4) placed relative to that point.
      const Point<3> corner_points[dim] = {Point<dim>((corners[1][0]-corners[0][0])*constants::degree_to_radians,
                                                      (corners[1][1]-corners[0][1])*constants::degree_to_radians,
                                                      0),
                                           Point<dim>((corners[3][0]-corners[0][0])*constants::degree_to_radians,
                                                      (corners[3][1]-corners[0][1])*constants::degree_to_radians,
                                                      0),
                                           Point<dim>(0,
                                                      0,
                                                      bottom_depth)
                                          };
      const unsigned int  subdivisions[dim] = {EW_subdiv, NS_subdiv, depth_subdiv};

      GridGenerator::subdivided_parallelepiped (coarse_grid, subdivisions, corner_points, true);

      // Shift the grid point at (0,0,0) (and the rest of the
      // points with it) to the correct location at corner[0] at a
      // negative depth.
      const Point<3> base_point(corners[0][0] * constants::degree_to_radians,
                                corners[0][1] * constants::degree_to_radians,
                                -bottom_depth);
      GridTools::shift(base_point,coarse_grid);

      // Transform to the ellipsoid surface
      GridTools::transform (
        [&](const Point<dim> &x) -> Point<dim>
      {
        return manifold->push_forward(x);
      },
      coarse_grid);

      // Also attach the real manifold to slot 15 (arbitrarily chosen).
      // We won't use it during regular operation, but we set manifold_ids for
      // all cells, faces and edges immediately before refinement and
      // clear it again afterwards
      coarse_grid.set_manifold (15, *manifold);

      coarse_grid.signals.pre_refinement.connect (
        [&] {set_manifold_ids(coarse_grid);});
      coarse_grid.signals.post_refinement.connect (
        [&] {clear_manifold_ids(coarse_grid);});
      coarse_grid.signals.post_refinement.connect (
        [&] {this->set_boundary_ids(coarse_grid);});
    }

    template <int dim>
    void
    EllipsoidalChunk<dim>::create_coarse_mesh(parallel::distributed::Triangulation<dim> &/*coarse_grid*/) const
    {
      Assert(false, ExcNotImplemented());
    }

    template <int dim>
    void
    EllipsoidalChunk<dim>::set_boundary_ids(parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      // set all boundary indicators. we want the edges to be curved
      // as well. for this, it matters in which order we call
      // set_all_boundary_indicators() -- we have to do it last for
      // the inner and outer boundary, which conveniently is what
      // happens in the following loop
      for (const auto &cell : coarse_grid.active_cell_iterators())
        for (const unsigned int f : cell->face_indices())
          if (cell->face(f)->at_boundary())
            cell->face(f)->set_boundary_id(f);
    }

    template <int dim>
    std::set<types::boundary_id>
    EllipsoidalChunk<dim>::get_used_boundary_indicators() const
    {
      // boundary indicators are zero through 2*dim-1
      std::set<types::boundary_id> s;
      for (unsigned int i = 0; i < 2 * dim; ++i)
        s.insert(i);
      return s;
    }

    template <int dim>
    std::map<std::string,types::boundary_id>
    EllipsoidalChunk<dim>::
    get_symbolic_boundary_names_map () const
    {
      switch (dim)
        {
          case 2:
          {
            return
            {
              {"east",   0},
              {"west",   1},
              {"bottom", 2},
              {"top",    3}
            };
          }

          case 3:
          {
            return
            {
              {"east",   0},
              {"west",   1},
              {"north",  2},
              {"south",  3},
              {"bottom", 4},
              {"top",    5}
            };
          }
        }

      Assert (false, ExcNotImplemented());
      return {};
    }

    template <int dim>
    void
    EllipsoidalChunk<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Ellipsoidal chunk");
        {
          prm.declare_entry("NE corner",
                            "",
                            Patterns::Anything(),
                            "Longitude:latitude in degrees of the North-East corner point of model region."
                            "The North-East direction is positive. If one of the three corners is not provided "
                            "the missing corner value will be calculated so all faces are parallel.");
          prm.declare_entry("NW corner",
                            "",
                            Patterns::Anything(),
                            "Longitude:latitude in degrees of the North-West corner point of model region. "
                            "The North-East direction is positive. If one of the three corners is not provided "
                            "the missing corner value will be calculated so all faces are parallel.");
          prm.declare_entry("SW corner",
                            "",
                            Patterns::Anything(),
                            "Longitude:latitude in degrees of the South-West corner point of model region. "
                            "The North-East direction is positive. If one of the three corners is not provided "
                            "the missing corner value will be calculated so all faces are parallel.");
          prm.declare_entry("SE corner",
                            "",
                            Patterns::Anything(),
                            "Longitude:latitude in degrees of the South-East corner point of model region. "
                            "The North-East direction is positive. If one of the three corners is not provided "
                            "the missing corner value will be calculated so all faces are parallel.");
          prm.declare_entry("Depth",
                            "500000.0",
                            Patterns::Double(0.),
                            "Bottom depth of model region.");
          prm.declare_entry("Semi-major axis",
                            "6378137.0",
                            Patterns::Double(0.),
                            "The semi-major axis (a) of an ellipsoid. This is the radius for a sphere (eccentricity=0). Default WGS84 semi-major axis.");
          prm.declare_entry("Eccentricity",
                            "8.1819190842622e-2",
                            Patterns::Double(0.),
                            "Eccentricity of the ellipsoid. Zero is a perfect sphere, default (8.1819190842622e-2) is WGS84.");
          prm.declare_entry("East-West subdivisions",
                            "1",
                            Patterns::Integer(0),
                            "The number of subdivisions of the coarse (initial) mesh in the East-West direction.");
          prm.declare_entry("North-South subdivisions",
                            "1",
                            Patterns::Integer(0),
                            "The number of subdivisions of the coarse (initial) mesh in the North-South direction.");
          prm.declare_entry("Depth subdivisions",
                            "1",
                            Patterns::Integer(0),
                            "The number of subdivisions of the coarse (initial) mesh in depth.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    EllipsoidalChunk<dim>::parse_parameters(ParameterHandler &prm)
    {
      AssertThrow (dim == 3, ExcMessage("This geometry can currently only be used in 3d."));

      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Ellipsoidal chunk");
        {
          // Get latitude and longitudes defining region of interest from
          // the parameter file.
          corners.resize(4);
          std::string NEcorner = prm.get("NE corner");
          std::string NWcorner = prm.get("NW corner");
          std::string SWcorner = prm.get("SW corner");
          std::string SEcorner = prm.get("SE corner");

          // make a list of what corners are present and check that there should be one or two corner missing,
          // otherwise throw an exception
          std::vector<bool> present(4,true);
          unsigned int missing = 0;
          if (NEcorner == "")
            {
              ++missing;
              present[0] = false;
            }
          if (NWcorner == "")
            {
              ++missing;
              present[1] = false;
            }
          if (SWcorner == "")
            {
              ++missing;
              present[2] = false;
            }
          if (SEcorner == "")
            {
              ++missing;
              present[3] = false;
            }

          AssertThrow (missing != 0,
                       ExcMessage ("Only two or three of the four corner points should be provided."));
          AssertThrow (missing == 1 || missing == 2,
                       ExcMessage ("Please provide two or three corner points."));

          std::vector<double> temp;

          if (present[0])
            {
              temp = Utilities::string_to_double(Utilities::split_string_list(NEcorner,':'));
              AssertThrow (temp.size() == 2, ExcMessage ("Two coordinates should be given for the NE-corner (longitude:latitude)."));
              corners[0] = Point<2>(temp[0],temp[1]);
            }
          else
            corners[0]= Point<2>();

          if (present[1])
            {
              temp = Utilities::string_to_double(Utilities::split_string_list(NWcorner,':'));
              AssertThrow (temp.size() == 2, ExcMessage ("Two coordinates should be given for the NW-corner (longitude:latitude)."));
              corners[1] = Point<2>(temp[0],temp[1]);
            }
          else
            corners[1]= Point<2>();

          if (present[2])
            {
              temp = Utilities::string_to_double(Utilities::split_string_list(SWcorner,':'));
              AssertThrow (temp.size() == 2, ExcMessage ("Two coordinates should be given for the SW-corner (longitude:latitude)."));
              corners[2] = Point<2>(temp[0],temp[1]);
            }
          else
            corners[2]= Point<2>();

          if (present[3])
            {
              temp = Utilities::string_to_double(Utilities::split_string_list(SEcorner,':'));
              AssertThrow (temp.size() == 2, ExcMessage ("Two coordinates should be given for the SE-corner (longitude:latitude)."));
              corners[3] = Point<2>(temp[0],temp[1]);
            }
          else
            corners[3]= Point<2>();


          bottom_depth = prm.get_double("Depth");
          semi_major_axis_a = prm.get_double("Semi-major axis");
          eccentricity = prm.get_double("Eccentricity");
          semi_minor_axis_b = std::sqrt((1 - Utilities::fixed_power<2>(eccentricity)) * Utilities::fixed_power<2>(semi_major_axis_a));
          EW_subdiv = prm.get_integer("East-West subdivisions");
          NS_subdiv = prm.get_integer("North-South subdivisions");
          depth_subdiv = prm.get_integer("Depth subdivisions");

          // Check whether the corners of the rectangle are really placed correctly
          if (present[0] == true && present[1] == true)
            {
              AssertThrow (corners[0][0] >= corners[1][0],
                           ExcMessage ("The longitude of the NE corner (" + boost::lexical_cast<std::string>(corners[0][0]) + ") cannot be smaller than the longitude of the NW corner (" + boost::lexical_cast<std::string>(corners[1][0]) + ")."));
            }
          if (present[0] == true && present[3] == true)
            {
              AssertThrow (corners[0][1] >= corners[3][1],
                           ExcMessage ("The latitude of the NE (" + boost::lexical_cast<std::string>(corners[0][1]) + ") corner cannot be smaller than the latitude of the SE corner (" + boost::lexical_cast<std::string>(corners[3][1]) + ")."));
            }
          if (present[2] == true && present[3] == true)
            {
              AssertThrow (corners[2][0] <= corners[3][0],
                           ExcMessage ("The longitude of the SW (" + boost::lexical_cast<std::string>(corners[2][0]) + ") corner cannot be larger than the longitude of the SE corner (" + boost::lexical_cast<std::string>(corners[3][0]) + ")."));
            }
          if (present[2] == true && present[1] == true)
            {
              AssertThrow (corners[2][1] <= corners[1][1],
                           ExcMessage ("The latitude of the SW corner (" + boost::lexical_cast<std::string>(corners[2][1]) + ") cannot be larger than the latitude of the NW corner (" + boost::lexical_cast<std::string>(corners[1][1]) + ")."));
            }
          if (missing == 2)
            {
              AssertThrow ((present[0] == true && present[2] == true) || (present[1] == true && present[3] == true),
                           ExcMessage ("Please provide two opposing corners."));

              if (present[0] == true && present[2] == true)
                AssertThrow (corners[0][0] > corners[2][0] && corners[0][1] > corners[2][1],
                             ExcMessage ("The Northeast corner (" + boost::lexical_cast<std::string>(corners[0][0]) + ":"  + boost::lexical_cast<std::string>(corners[0][1]) + ") must be strictly north and east of the Southwest corner (" + boost::lexical_cast<std::string>(corners[2][0]) + ":"  + boost::lexical_cast<std::string>(corners[2][1]) + ") when only two points are given."));

              if (present[1] == true && present[3] == true)
                AssertThrow (corners[1][0] < corners[3][0] && corners[1][1] > corners[3][1],
                             ExcMessage ("The Northwest corner (" + boost::lexical_cast<std::string>(corners[1][0]) + ":"  + boost::lexical_cast<std::string>(corners[1][1]) + ") must be strictly north and west of the Southeast corner (" + boost::lexical_cast<std::string>(corners[3][0]) + ":"  + boost::lexical_cast<std::string>(corners[3][1]) + ") when only two points are given."));
            }


          // If one or two of the corners is not provided, calculate it.
          if (missing == 1)
            {
              for (unsigned int i = 0; i <= 3; ++i)
                {
                  // find the corner which is missing, if there is one
                  if (present[i] == false)
                    {
                      // find the location of this corner
                      // the values for the values is:
                      // [i] = [i+1] + ([i+3] - [i+2])
                      unsigned int ip = i + 1;
                      unsigned int ipp = i + 2;
                      unsigned int ippp = i + 3;

                      if (ip > 3)
                        ip -= 4;
                      if (ipp > 3)
                        ipp -= 4;
                      if (ippp > 3)
                        ippp -=4;

                      corners[i] = corners[ip];
                      corners[i][0] = corners[ip][0] + (corners[ippp][0] - corners[ipp][0]);
                      corners[i][1] = corners[ip][1] + (corners[ippp][1] - corners[ipp][1]);
                    }
                }
            }
          else if (missing == 2)
            {
              // find the location of this corner
              if (present[0] == false)
                {
                  corners[0] = corners[1];
                  corners[2] = corners[3];
                  corners[0][0] = corners[3][0];
                  corners[2][0] = corners[1][0];
                }
              else
                {
                  corners[1] = corners[0];
                  corners[3] = corners[2];
                  corners[1][0] = corners[2][0];
                  corners[3][0] = corners[0][0];
                }
            }
          // Check that the calculated corners also obey the rules for the location of the corners.
          Assert (corners[0][0] >= corners[1][0],
                  ExcMessage ("The longitude of the NE corner (" + boost::lexical_cast<std::string>(corners[0][0]) + ") cannot be smaller than the longitude of the NW corner (" + boost::lexical_cast<std::string>(corners[1][0]) + "). This is an internal check, if you see this please contact the developer."));

          Assert (corners[0][1] >= corners[3][1],
                  ExcMessage ("The latitude of the NE (" + boost::lexical_cast<std::string>(corners[0][1]) + ") corner cannot be smaller than the latitude of the SE corner (" + boost::lexical_cast<std::string>(corners[3][1]) + "). This is an internal check, if you see this please contact the developer."));

          Assert (corners[2][0] <= corners[3][0],
                  ExcMessage ("The longitude of the SW (" + boost::lexical_cast<std::string>(corners[2][0]) + ") corner cannot be larger than the longitude of the SE corner (" + boost::lexical_cast<std::string>(corners[3][0]) + "). This is an internal check, if you see this please contact the developer."));

          Assert (corners[2][1] <= corners[1][1],
                  ExcMessage ("The latitude of the SW corner (" + boost::lexical_cast<std::string>(corners[2][1]) + ") cannot be larger than the latitude of the NW corner (" + boost::lexical_cast<std::string>(corners[1][1]) + "). This is an internal check, if you see this please contact the developer."));

          westLongitude = corners[2][0];
          eastLongitude = corners[0][0];
          northLatitude = corners[0][1];
          southLatitude = corners[2][1];
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    double
    EllipsoidalChunk<dim>::depth(const Point<dim> &position) const
    {
      return std::max(std::min(-manifold->pull_back(position)[dim-1], maximal_depth()), 0.0);
    }

    template <int dim>
    double
    EllipsoidalChunk<dim>::height_above_reference_surface(const Point<dim> &/*position*/) const
    {
      AssertThrow(false, ExcMessage("Function height_above_reference_surface is not yet implemented "
                                    "for the ellipsoidal chunk geometry model. "
                                    "Consider using a box, spherical shell, or chunk.") );
      return numbers::signaling_nan<double>();
    }


    template <int dim>
    double
    EllipsoidalChunk<dim>::maximal_depth() const
    {
      return bottom_depth;
    }

    template <int dim>
    double
    EllipsoidalChunk<dim>::get_radius(const Point<dim> &position) const
    {
      const Point<dim> long_lat_depth = manifold->pull_back(position);
      return semi_major_axis_a / (std::sqrt(1 - eccentricity * eccentricity * std::sin(long_lat_depth[1]) * std::sin(long_lat_depth[1])));
    }

    template <int dim>
    double
    EllipsoidalChunk<dim>::get_eccentricity() const
    {
      return eccentricity;
    }

    template <int dim>
    const std::vector<Point<2>> &
    EllipsoidalChunk<dim>::get_corners() const
    {
      return corners;
    }

    template <int dim>
    double
    EllipsoidalChunk<dim>::get_semi_minor_axis_b() const
    {
      return semi_minor_axis_b;
    }

    template <int dim>
    double
    EllipsoidalChunk<dim>::get_semi_major_axis_a() const
    {
      return semi_major_axis_a;
    }


    template <int dim>
    double
    EllipsoidalChunk<dim>::length_scale() const
    {
      // As described in the first ASPECT paper, a length scale of
      // 10km = 1e4m works well for the pressure scaling for earth
      // sized spherical shells. use a length scale that
      // yields this value for the R0,R1 corresponding to earth
      // but otherwise scales like (R1-R0)
      return 1e4 * maximal_depth() / (6336000.-3481000.);
    }

    template <int dim>
    Point<dim>
    EllipsoidalChunk<dim>::representative_point(const double /*depth*/) const
    {
      return Point<dim>();
    }

    template <>
    Point<3>
    EllipsoidalChunk<3>::representative_point(const double depth) const
    {
      const int dim = 3;
      Assert(depth >= 0, ExcMessage("Given depth must be positive or zero."));
      Assert(depth <= maximal_depth(),
             ExcMessage("Given depth must be less than or equal to the maximal depth of this geometry."));

      // Choose a point on the center axis of the domain
      Point<dim> p = Point<3>((eastLongitude + westLongitude) * 0.5 * constants::degree_to_radians,
                              (southLatitude + northLatitude) * 0.5 * constants::degree_to_radians,
                              -depth);

      return manifold->push_forward(p);
    }


    template <int dim>
    bool
    EllipsoidalChunk<dim>::point_is_in_domain(const Point<dim> &point) const
    {
      AssertThrow(!this->get_parameters().mesh_deformation_enabled ||
                  this->get_timestep_number() == 0,
                  ExcMessage("After displacement of the mesh, this function can no longer be used to determine whether a point lies in the domain or not."));

      // dim = 3
      const Point<dim> ellipsoidal_point = manifold->pull_back(point);

      // compare deflection from the ellipsoid surface
      if (ellipsoidal_point[dim-1] > 0.0+std::numeric_limits<double>::epsilon()*bottom_depth ||
          -ellipsoidal_point[dim-1] > bottom_depth+std::numeric_limits<double>::epsilon()*bottom_depth)
        return false;

      // compare lon/lat
      if (!Utilities::polygon_contains_point<dim>(corners, Point<2>(ellipsoidal_point[0]*constants::radians_to_degree,ellipsoidal_point[1]*constants::radians_to_degree)))
        return false;

      return true;
    }

    template <int dim>
    std::array<double,dim>
    EllipsoidalChunk<dim>::cartesian_to_natural_coordinates(const Point<dim> &position_point) const
    {
      // the chunk manifold works internally with a vector with longitude, latitude, depth.
      // We need to output radius, longitude, latitude to be consistent.
      // Ignore the topography by calling pull_back_ellipsoid to avoid a loop when calling the
      // AsciiDataBoundary for topography which uses this function....
      Point<3> cartesian_point;
      for (unsigned int d=0; d<dim; ++d)
        cartesian_point[d] = position_point[d];

      Point<3> transformed_point = manifold->pull_back_ellipsoid(cartesian_point, semi_major_axis_a, eccentricity);

      const double radius =  semi_major_axis_a /
                             (std::sqrt(1 - eccentricity * eccentricity * std::sin(transformed_point[1]) * std::sin(transformed_point[1])));
      std::array<double,dim> position_array;
      position_array[0] = radius + transformed_point(2);
      position_array[1] = transformed_point(0);
      position_array[2] = transformed_point(1);

      return position_array;
    }


    template <int dim>
    aspect::Utilities::Coordinates::CoordinateSystem
    EllipsoidalChunk<dim>::natural_coordinate_system() const
    {
      return aspect::Utilities::Coordinates::CoordinateSystem::ellipsoidal;
    }


    template <>
    Point<3>
    EllipsoidalChunk<3>::natural_to_cartesian_coordinates(const std::array<double,3> &position_tensor) const
    {
      // We receive radius, longitude, latitude and we need to turn it first back into
      // longitude, latitude, depth for internal use, and push_forward to cartesian coordinates.
      // Ignore the topography by calling push_forward_ellipsoid to avoid a loop when calling the
      // AsciiDataBoundary for topography which uses this function....
      Point<3> position_point;
      position_point(0) = position_tensor[1];
      position_point(1) = position_tensor[2];

      const double radius = semi_major_axis_a / (std::sqrt(1 - eccentricity * eccentricity * std::sin(position_point(1)) * std::sin(position_point(1))));
      position_point(2) = position_tensor[0] - radius;

      Point<3> transformed_point = manifold->push_forward_ellipsoid(position_point, semi_major_axis_a, eccentricity);
      return transformed_point;
    }


    template <>
    Point<2>
    EllipsoidalChunk<2>::natural_to_cartesian_coordinates(const std::array<double,2> &/*position_tensor*/) const
    {
      Assert(false, ExcMessage("This geometry model doesn't support 2d."));
      return {};
    }

  }
}

template <int dim>
typename aspect::GeometryModel::internal::EllipsoidalChunkGeometry<dim>
aspect::GeometryModel::EllipsoidalChunk<dim>::get_manifold() const
{
  return *manifold;
}

// explicit instantiations
namespace aspect
{
  namespace GeometryModel
  {
    namespace internal
    {
#define INSTANTIATE(dim) \
  template class EllipsoidalChunkGeometry<dim>;
      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }

    ASPECT_REGISTER_GEOMETRY_MODEL(EllipsoidalChunk,
                                   "ellipsoidal chunk",
                                   "A 3d chunk geometry that accounts for Earth's ellipticity (default assuming the "
                                   "WGS84 ellipsoid definition) which can be defined in non-coordinate directions. "
                                   "In the description of the ellipsoidal chunk, two of the ellipsoidal axes have the "
                                   "same length so that there is only a semi-major axis and a semi-minor axis. "
                                   "The user has two options for creating an ellipsoidal chunk geometry: 1) by defining "
                                   "two opposing points (SW and NE or NW and SE) a coordinate parallel ellipsoidal "
                                   "chunk geometry will be created. 2) by defining three points a non-coordinate "
                                   "parallel ellipsoidal chunk will be created. The points are defined in the input "
                                   "file by longitude:latitude. It is also possible to define additional subdivisions "
                                   "of the mesh in each direction. The boundary of the domain is formed by linear "
                                   "interpolation in longitude-latitude space between adjacent points "
                                   "(i.e. $[lon, lat](f) = [lon1 \\cdot f + lon2 \\cdot(1-f), lat1 \\cdot f + lat2 \\cdot (1-f)]$, "
                                   "where f is a value between 0 and 1). Faces of the model are defined as "
                                   "0, west; 1,east; 2, south; 3, north; 4, inner; 5, outer.\n\n"
                                   "This geometry model supports initial topography for deforming the initial mesh.")
  }
}
