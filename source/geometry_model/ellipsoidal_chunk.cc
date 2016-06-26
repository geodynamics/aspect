/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/geometry_model/ellipsoidal_chunk.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <boost/lexical_cast.hpp>
#include <aspect/compat.h>



/**
 * This geometry model implements an (3d) ellipsoidal_chunk geometry which can be non-coordinate parallel.
 * @author This plugin is a joined effort of Menno Fraters, D Sarah Stamps and Wolfgang Bangerth
 */

namespace aspect
{
  namespace GeometryModel
  {
    using namespace dealii;

    template <int dim>
    bool
    EllipsoidalChunk<dim>::In2dPolygon(dealii::Point<2> &point,const std::vector<std::vector<double> > &pointList)
    {
      /**
       * This code has been based on http://geomalgorithms.com/a03-_inclusion.html,
       * and therefore requires the following copyright notice:
       *
       * Copyright 2000 softSurfer, 2012 Dan Sunday
       * This code may be freely used and modified for any purpose
       * providing that this copyright notice is included with it.
       * SoftSurfer makes no warranty for this code, and cannot be held
       * liable for any real or imagined damage resulting from its use.
       * Users of this code must verify correctness for their application.
       */
      int pointNo = pointList.size();
      int    wn = 0;    // the  winding number counter
      int   j=pointNo-1;


      // loop through all edges of the polygon
      for (int i=0; i<pointNo; i++)
        {
          // edge from V[i] to  V[i+1]
          if (pointList[j][1] <= point[1])
            {
              // start y <= P.y
              if (pointList[i][1]  > point[1])      // an upward crossing
                if (( (pointList[i][0] - pointList[j][0]) * (point[1] - pointList[j][1])
                      - (point[0] -  pointList[j][0]) * (pointList[i][1] - pointList[j][1]) ) > 0)
                  {
                    // P left of  edge
                    ++wn;            // have  a valid up intersect
                  }
            }
          else
            {
              // start y > P.y (no test needed)
              if (pointList[i][1]  <= point[1])     // a downward crossing
                if (( (pointList[i][0] - pointList[j][0]) * (point[1] - pointList[j][1])
                      - (point[0] -  pointList[j][0]) * (pointList[i][1] - pointList[j][1]) ) < 0)
                  {
                    // P right of  edge
                    --wn;            // have  a valid down intersect
                  }
            }
          j=i;
        }
      return (wn != 0);
    }

    /**
     * the EllipsoidalChunkTopography class
     */
    // constructor
    template <int dim>
    EllipsoidalChunk<dim>::EllipsoidalChunkTopography::EllipsoidalChunkTopography()
    {}


    template <int dim>
    double
    EllipsoidalChunk<dim>::EllipsoidalChunkTopography::value (const double lon,
                                                              const double lat) const
    {
      double value = 0;

      switch (topo_type)
        {
          case NO_TOPOGRAPHY:
            return 0;
            break;

          case PRM_EXACT:
            for (unsigned int i = 0; i < point_lists.size(); i++)
              {
                Point<2> p1(lon * 180/numbers::PI,lat * 180/numbers::PI);

                if (In2dPolygon(p1, point_lists[i]))
                  {
                    value = topography_values[i];
                  }
              }
            return value;
            break;

          case PRM_UNIFORM_GRID_INTERPOLATED:
            if (topography_data != NULL)
              return static_cast<Functions::InterpolatedUniformGridData<2>*>(topography_data)->value (Point<2>(lat * 180/numbers::PI,
                     lon * 180/numbers::PI));
            break;

          case FILE_UNIFORM_GRID:
            return static_cast<Functions::InterpolatedUniformGridData<2>*>(topography_data)->value (Point<2>(lat * 180/numbers::PI,
                   lon * 180/numbers::PI));
            break;

          case FILE_NONUNIFORM_GRID:
            return static_cast<Functions::InterpolatedTensorProductGridData<2>*>(topography_data)->value (Point<2>(lat * 180/numbers::PI,
                   lon * 180/numbers::PI));
            break;

          default:
            AssertThrow(false,ExcMessage ("This topography function for enum with value " + boost::lexical_cast<std::string>(topo_type) + " has not been implemented."));
            return 0;
            break;
        }
      Assert(false,ExcMessage ("The code should never reach this part. Please check the code or contact the developer."));
      return 0;
    }

    template <int dim>
    void
    EllipsoidalChunk<dim>::EllipsoidalChunkTopography::set_corners(std::vector<Point<2> > set_corners)
    {
      corners = set_corners;
    }

    template <int dim>
    void
    EllipsoidalChunk<dim>::EllipsoidalChunkTopography::set_topography_type(TopoTypes set_topo_type)
    {
      topo_type = set_topo_type;
    }

    template <int dim>
    void
    EllipsoidalChunk<dim>::EllipsoidalChunkTopography::set_point_lists(std::vector<std::vector<std::vector<double> > > set_point_lists)
    {
      point_lists = set_point_lists;
    }

    template <int dim>
    void
    EllipsoidalChunk<dim>::EllipsoidalChunkTopography::set_topography_values(std::vector<double> set_topography_values)
    {
      topography_values = set_topography_values;
    }


    template <int dim>
    void
    EllipsoidalChunk<dim>::EllipsoidalChunkTopography::set_topography_data (Function<2> *set_topography_data)
    {
      topography_data = set_topography_data;
    }


    template <int dim>
    void
    EllipsoidalChunk<dim>::EllipsoidalChunkTopography::set_grid_number_data_points(std::vector<double> &set_grid_number_data_points)
    {
      grid_number_data_points = set_grid_number_data_points;
    }

    template <int dim>
    void
    EllipsoidalChunk<dim>::EllipsoidalChunkTopography::set_topography_file(std::string set_topo_file_name)
    {
      topo_file = set_topo_file_name;
    }

    template <int dim>
    TopoTypes
    EllipsoidalChunk<dim>::EllipsoidalChunkTopography::get_topo_type () const
    {
      return topo_type;
    }

    template <int dim>
    std_cxx11::array<std::pair<double,double>,2>
    EllipsoidalChunk<dim>::EllipsoidalChunkTopography::get_endpoints () const
    {
      std_cxx11::array<std::pair<double,double>,2> endpoints;
      // Minimal and maximal longitude
      endpoints[0] = std::make_pair (corners[1][0]*numbers::PI/180, corners[0][0]*numbers::PI/180);
      // Minimal and maximal latitude
      endpoints[1] = std::make_pair (corners[3][1]*numbers::PI/180, corners[0][1]*numbers::PI/180);

#ifdef DEBUG
      for (unsigned int i=0; i<2; i++)
        Assert (endpoints[i].first < endpoints[i].second,
                ExcMessage ("The interval in each coordinate direction needs "
                            "to have positive size"));
#endif

      return endpoints;
    }

    template <int dim>
    std_cxx11::array<unsigned int,2>
    EllipsoidalChunk<dim>::EllipsoidalChunkTopography::get_number_of_intervals() const
    {
      std_cxx11::array<unsigned int,2> interval;
      interval[0] = grid_number_data_points[0]-1;
      interval[1] = grid_number_data_points[1]-1;
      return interval;
    }

    template <int dim>
    std::vector<double>
    EllipsoidalChunk<dim>::EllipsoidalChunkTopography::get_data () const
    {
      if (topo_type == PRM_UNIFORM_GRID_INTERPOLATED)
        {
          // TODO: use grid_data for this?
          // Fill grid based on polygons from prm
          std::vector<double> data(grid_number_data_points[0] * grid_number_data_points[1],0);
          double d_long = (corners[3][0]-corners[2][0])/grid_number_data_points[0];
          double d_lat = (corners[1][1]-corners[2][1])/grid_number_data_points[1];

          for (unsigned int i_long = 0; i_long < grid_number_data_points[0]; i_long++)
            {
              for (unsigned int i_lat = 0; i_lat < grid_number_data_points[1]; i_lat++)
                {
                  for (unsigned int i = 0; i < point_lists.size(); i++)
                    {
                      Point<2> p1((i_long * d_long + corners[2][0]) * 180/numbers::PI,(i_lat * d_lat + corners[2][1]) * 180/numbers::PI);

                      if (In2dPolygon(p1, point_lists[i]))
                        {
                          data[i_long * grid_number_data_points[1] + i_lat] = topography_values[i];
                        }
                    }
                }
            }
          return data;
        }
      else if (topo_type == FILE_UNIFORM_GRID || topo_type == FILE_NONUNIFORM_GRID)
        {
          // Return the data from file
          return grid_data;
        }
      else
        {
          AssertThrow(false,ExcNotImplemented());

        }
      return std::vector<double>();



    }

    template <int dim>
    void
    EllipsoidalChunk<dim>::EllipsoidalChunkTopography::get_data_from_file ()
    {
      unsigned int n_data_points = grid_number_data_points[0] * grid_number_data_points[1];

      grid_data.resize(n_data_points,0);
      std::vector<double> lon_points;
      std::vector<double> lat_points;

      // In file stream
      std::ifstream in_topo(topo_file.c_str(), std::ios::in);
      // Check whether file exists, if not, throw exception
      AssertThrow (in_topo,
                   ExcMessage (std::string("Couldn't open file ") + topo_file));

      double topo=0.0;

      // Read in coordinates of nonuniform grid
      // TODO: format of nonuniform grid file might change
      if (topo_type == FILE_NONUNIFORM_GRID)
        {
          // First longitude coordinates
          for (unsigned int i=0; i<grid_number_data_points[0]; i++)
            {
              if (!(in_topo >> topo))
                {
                  AssertThrow(false, ExcMessage("Could not read longitude point " + dealii::Utilities::int_to_string(i)
                                                + " of file " + topo_file));
                }
              if (i >= 1)
                AssertThrow(topo > lon_points[i-1], ExcMessage("Longitude grid coordinates must be strictly ascending."));
              lon_points.push_back(topo);
            }
          // Then latitude coordinates
          for (unsigned int i=0; i<grid_number_data_points[1]; i++)
            {
              if (!(in_topo >> topo))
                {
                  AssertThrow(false, ExcMessage("Could not read latitude point " + dealii::Utilities::int_to_string(i)
                                                + " of file " + topo_file));
                }
              if (i >= 1)
                AssertThrow(topo > lat_points[i-1], ExcMessage("Latitude grid coordinates must be strictly ascending."));
              lat_points.push_back(topo);
            }

          nonuniform_grid_coordinates[0]=lon_points;
          nonuniform_grid_coordinates[1]=lat_points;
        }


      // Read in topography values of uniform or nonuniform grid
      for (unsigned int i=0; i<n_data_points; i++)
        {
          if (!(in_topo >> topo))
            {
              AssertThrow(false, ExcMessage("Could not read point " + dealii::Utilities::int_to_string(i)
                                            + " of file " + topo_file));

            }
          grid_data[i]=topo;
        }

    }


    template <int dim>
    std_cxx11::array< std::vector< double >, 2 >
    EllipsoidalChunk<dim>::EllipsoidalChunkTopography::get_coordinates() const
    {
      Assert (topo_type == FILE_NONUNIFORM_GRID, ExcMessage("Nonuniform grid not selected, so no need for its coordinates. "));
      return nonuniform_grid_coordinates;
    }


    /**
     * the EllipsoidalChunkGeometry class
     */

// constructor
    template <int dim>
    EllipsoidalChunk<dim>::EllipsoidalChunkGeometry::EllipsoidalChunkGeometry()
      :
      semi_major_axis_a (-1),
      eccentricity (-1),
      semi_minor_axis_b (-1),
      bottom_depth (-1)
    {}

    template <int dim>
    void
    EllipsoidalChunk<dim>::EllipsoidalChunkGeometry::set_manifold_parameters(const double para_semi_major_axis_a,
                                                                             const double para_eccentricity,
                                                                             const double para_semi_minor_axis_b,
                                                                             const double para_bottom_depth,
                                                                             const std::vector<Point<2> > &para_corners)
    {
      semi_major_axis_a = para_semi_major_axis_a;
      eccentricity = para_eccentricity;
      semi_minor_axis_b = para_semi_minor_axis_b;
      bottom_depth = para_bottom_depth;
      corners=para_corners;
    }

    template <int dim>
    Point<3>
    EllipsoidalChunk<dim>::EllipsoidalChunkGeometry::push_forward_ellipsoid(const Point<3> &phi_theta_d, const double semi_major_axis_a, const double eccentricity) const
    {
      const double phi   = phi_theta_d[0];
      const double theta = phi_theta_d[1];
      const double d     = phi_theta_d[2];

      const double R_bar = semi_major_axis_a / std::sqrt(1 - (eccentricity * eccentricity *
                                                              std::sin(theta) * std::sin(theta)));

      return Point<3> ((R_bar + d) * std::cos(phi) * std::cos(theta),
                       (R_bar + d) * std::sin(phi) * std::cos(theta),
                       ((1 - eccentricity * eccentricity) * R_bar + d) * std::sin(theta));
    }

    template <int dim>
    Point<3>
    EllipsoidalChunk<dim>::EllipsoidalChunkGeometry::pull_back_ellipsoid(const Point<3> &x, const double semi_major_axis_a, const double eccentricity) const
    {
      const double R    = semi_major_axis_a;
      const double b      = std::sqrt(R * R * (1 - eccentricity * eccentricity));
      const double ep     = std::sqrt((R * R - b * b) / (b * b));
      const double p      = std::sqrt(x(0) * x(0) + x(1) * x(1));
      const double th     = std::atan2(R * x(2), b * p);
      const double phi    = std::atan2(x(1), x(0));
      const double theta  = std::atan2(x(2) + ep * ep * b * std::pow(std::sin(th),3),
                                       (p - (eccentricity * eccentricity * R  * std::pow(std::cos(th),3))));
      const double R_bar = R / (std::sqrt(1 - eccentricity * eccentricity * std::sin(theta) * std::sin(theta)));
      const double R_plus_d = p / std::cos(theta);

      Point<3> phi_theta_d;
      if (phi < -numbers::PI)
        phi_theta_d[0] = phi + 2*numbers::PI;
      else if (phi > numbers::PI)
        phi_theta_d[0] = phi - 2*numbers::PI;
      else
        phi_theta_d[0] = phi;

      phi_theta_d[1] = theta;
      phi_theta_d[2] = R_plus_d - R_bar;
      return phi_theta_d;
    }

    // From topo depth to normal depth
    template <int dim>
    Point<3>
    EllipsoidalChunk<dim>::EllipsoidalChunkGeometry::push_forward_topography(const Point<3> &phi_theta_d_hat) const
    {
      const double d_hat = phi_theta_d_hat[2]; // long, lat, depth
      const double h     = topography.value(phi_theta_d_hat[0],
                                            phi_theta_d_hat[1]);

      const double d = d_hat + (d_hat + bottom_depth)/bottom_depth*h;
      const Point<3> phi_theta_d (phi_theta_d_hat[0],
                                  phi_theta_d_hat[1],
                                  d);
      return phi_theta_d;
    }

    // from normal depth to topo depth
    template <int dim>
    Point<3>
    EllipsoidalChunk<dim>::EllipsoidalChunkGeometry::pull_back_topography(const Point<3> &phi_theta_d) const
    {
      const double d = phi_theta_d[2];
      const double h = topography.value(phi_theta_d[0],
                                        phi_theta_d[1]);
      const double d_hat = bottom_depth * (d-h)/(bottom_depth+h);
      const Point<3> phi_theta_d_hat (phi_theta_d[0],
                                      phi_theta_d[1],
                                      d_hat);
      return phi_theta_d_hat;
    }

    /**
     * TODO: These functions (pull back and push forward) should be changed that they always
     * take and return 3D points, because 2D points make no sense for an ellipsoid, even with
     * a 2D triangulation. To do this correctly we need to add the spacedim to the triangulation
     * in ASPECT. What is now presented is just a temporary fix to get acces to the pull back
     * function from outside. The push forward function can't be fixed in this way, because
     * it is used by a bind statement.
     */
    template <int dim>
    Point<3>
    EllipsoidalChunk<dim>::EllipsoidalChunkGeometry::pull_back(const Point<3> &space_point) const
    {
      AssertThrow (dim == 3,ExcMessage ("This can not be done with 2D points."));
      return pull_back_topography(pull_back_ellipsoid (space_point, semi_major_axis_a, eccentricity));

    }

    template <int dim>
    Point<2>
    EllipsoidalChunk<dim>::EllipsoidalChunkGeometry::pull_back(const Point<2> &space_point) const
    {
      AssertThrow (dim == 3,ExcMessage ("This can not be done with 2D points.")); //TODO: check if space point is 2 or not.
      return space_point;

    }

    template <int dim>
    Point<3>
    EllipsoidalChunk<dim>::EllipsoidalChunkGeometry::push_forward(const Point<3> &chart_point) const
    {
      AssertThrow (dim == 3,ExcMessage ("This can not be done with 2D points."));

      return push_forward_ellipsoid (push_forward_topography(chart_point), semi_major_axis_a, eccentricity);
    }


    template <int dim>
    void
    EllipsoidalChunk<dim>::initialize ()
    {
      switch (manifold.topography.get_topo_type())
        {
          case NO_TOPOGRAPHY:
            break;

          case PRM_EXACT:
            break;

          case PRM_UNIFORM_GRID_INTERPOLATED:
            // create the uniform grid and add it to the pointer.
            manifold.topography.set_topography_data (new Functions::InterpolatedUniformGridData<2> (manifold.topography.get_endpoints(),
                                                     manifold.topography.get_number_of_intervals(),
                                                     Table<2,double> (manifold.topography.get_number_of_intervals()[0]+1,
                                                                      manifold.topography.get_number_of_intervals()[1]+1,
                                                                      manifold.topography.get_data().begin())));

            // TODO: add delete in destructor?
            break;
          case FILE_UNIFORM_GRID:
            // read in the uniform grid topography values and
            // add them to the pointer
            manifold.topography.get_data_from_file();
            manifold.topography.set_topography_data (new Functions::InterpolatedUniformGridData<2> (manifold.topography.get_endpoints(),
                                                     manifold.topography.get_number_of_intervals(),
                                                     Table<2,double> (manifold.topography.get_number_of_intervals()[0]+1,
                                                                      manifold.topography.get_number_of_intervals()[1]+1,
                                                                      manifold.topography.get_data().begin())));
            break;

          case FILE_NONUNIFORM_GRID:
          {
            // read in the nonuniform grid coordinates and topography values and
            // add them to the pointer
            manifold.topography.get_data_from_file();
            manifold.topography.set_topography_data (new Functions::InterpolatedTensorProductGridData<2> (manifold.topography.get_coordinates(),
                                                     Table<2,double> (manifold.topography.get_number_of_intervals()[0]+1,
                                                                      manifold.topography.get_number_of_intervals()[1]+1,
                                                                      manifold.topography.get_data().begin())));
          }
          break;

          default:
            AssertThrow(false,ExcMessage ("This topography function for enum with value " + boost::lexical_cast<std::string>(manifold.topography.get_topo_type()) + " has not been implemented."));
            break;
        }
    }

    template <>
    void
    EllipsoidalChunk<3>::create_coarse_mesh(parallel::distributed::Triangulation<3> &coarse_grid) const
    {
      const int dim = 3;

      /**
       * Generate parallelepiped grid with one point (point 0) at (0,0,0) and the
       * other corners (respectively corner 1, 2 and 3) placed relative to that point.
       */
      const Point<3> corner_points[dim] = {Point<dim>((corners[1][0]-corners[0][0])*numbers::PI/180,
                                                      (corners[1][1]-corners[0][1])*numbers::PI/180,
                                                      0),
                                           Point<dim>((corners[3][0]-corners[0][0])*numbers::PI/180,
                                                      (corners[3][1]-corners[0][1])*numbers::PI/180,
                                                      0),
                                           Point<dim>(0,
                                                      0,
                                                      bottom_depth)
                                          };

      const unsigned int  subdivisions[dim] = {EW_subdiv,NS_subdiv,depth_subdiv};

      GridGenerator::subdivided_parallelepiped (coarse_grid, subdivisions,corner_points, true);

      /**
       * Shift the grid point at (0,0,0) (and the rest of the
       * points with it) to the correct location at corner[0] at a
       * negative depth.
       */
      const Point<3> base_point(corners[0][0] *numbers::PI/180,corners[0][1] *numbers::PI/180,-bottom_depth);
      GridTools::shift(base_point,coarse_grid);

      // Transform to the ellipsoid surface
      GridTools::transform (std_cxx11::bind(&EllipsoidalChunk<3>::EllipsoidalChunkGeometry::push_forward,
                                            std_cxx11::cref(manifold),
                                            std_cxx11::_1),
                            coarse_grid);

      // also attach the real manifold to slot 15. we won't use it
      // during regular operation, but we set manifold_ids for all
      // cells, faces and edges immediately before refinement and
      // clear it again afterwards
      coarse_grid.set_manifold (15, manifold);

      coarse_grid.signals.pre_refinement.connect (std_cxx11::bind (&set_manifold_ids,
                                                                   std_cxx11::ref(coarse_grid)));
      coarse_grid.signals.post_refinement.connect (std_cxx11::bind (&clear_manifold_ids,
                                                                    std_cxx11::ref(coarse_grid)));
      coarse_grid.signals.post_refinement.connect(std_cxx11::bind (&EllipsoidalChunk<dim>::set_boundary_ids,
                                                                   std_cxx11::cref(*this),
                                                                   std_cxx11::ref(coarse_grid)));
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
      for (typename Triangulation<dim>::active_cell_iterator cell =
             coarse_grid.begin_active(); cell != coarse_grid.end(); ++cell)
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
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
            static const std::pair<std::string,types::boundary_id> mapping[]
              = { std::pair<std::string,types::boundary_id>("east",  0),
                  std::pair<std::string,types::boundary_id>("west",  1),
                  std::pair<std::string,types::boundary_id>("inner", 2),
                  std::pair<std::string,types::boundary_id>("outer", 3)
                };

            return std::map<std::string,types::boundary_id> (&mapping[0],
                                                             &mapping[sizeof(mapping)/sizeof(mapping[0])]);
          }

          case 3:
          {
            static const std::pair<std::string,types::boundary_id> mapping[]
              = { std::pair<std::string,types::boundary_id>("east",  0),
                  std::pair<std::string,types::boundary_id>("west",  1),
                  std::pair<std::string,types::boundary_id>("north", 2),
                  std::pair<std::string,types::boundary_id>("south", 3),
                  std::pair<std::string,types::boundary_id>("inner", 4),
                  std::pair<std::string,types::boundary_id>("outer", 5)
                };

            return std::map<std::string,types::boundary_id> (&mapping[0],
                                                             &mapping[sizeof(mapping)/sizeof(mapping[0])]);
          }
        }

      Assert (false, ExcNotImplemented());
      return std::map<std::string,types::boundary_id>();
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
                            "The North-East direction is positive. If one of the three corners is not provided"
                            "the missing corner value will be calculated so all faces are parallel.");
          prm.declare_entry("NW corner",
                            "",
                            Patterns::Anything(),
                            "Longitude:latitude in degrees of the North-West corner point of model region. "
                            "The North-East direction is positive. If one of the three corners is not provided"
                            "the missing corner value will be calculated so all faces are parallel.");
          prm.declare_entry("SW corner",
                            "",
                            Patterns::Anything(),
                            "Longitude:latitude in degrees of the South-West corner point of model region. "
                            "The North-East direction is positive. If one of the three corners is not provided"
                            "the missing corner value will be calculated so all faces are parallel.");
          prm.declare_entry("SE corner",
                            "",
                            Patterns::Anything(),
                            "Longitude:latitude in degrees of the South-East corner point of model region. "
                            "The North-East direction is positive. If one of the three corners is not provided"
                            "the missing corner value will be calculated so all faces are parallel.");
          prm.declare_entry("Depth",
                            "500000.0",
                            Patterns::Double(0),
                            "Bottom depth of model region.");
          prm.declare_entry("Semi-major axis",
                            "6378137.0",
                            Patterns::Double(0),
                            "The semi-major axis (a) of an ellipsoid. This is the radius for a sphere (eccentricity=0). Default WGS84 semi-major axis.");
          prm.declare_entry("Eccentricity",
                            "8.1819190842622e-2",
                            Patterns::Double(0),
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
          prm.enter_subsection("Topography");
          {
            prm.declare_entry("Topography type",
                              "No topography",
                              Patterns::Selection("No topography|From prm exact|From prm interpolated on uniform grid|From file uniform grid|From file nonuniform grid"),
                              "Select type of topography. "
                              "\n'No Topography' is the default option and creates no topography. "
                              "\n'From prm exact' uses the topography parameters directly to determine what elevation a certain point has. See the Topography parameters "
                              "for more information on how to define the topography function. "
                              "\n'From prm interpolated on uniform grid' uses the same function as 'From prm exact', but interpolates this information on a uniform grid "
                              "given by the parameters ."
                              "\n'From file uniform grid' reads in topography values defined on a uniform (in each direction) long,lat grid from file. "
                              "\n'From file nonuniform grid' reads in topography values defined on a variably spaced long,lat grid from file. " );
            prm.declare_entry("Topography parameters",
                              "",
                              Patterns::Anything(),
                              "Set the topography height, end with a |, and set the areas described by the points, separated by commas and coordinates separated by a ':'. "
                              "Seperate each topography feature by a semicolon. For example for two triangular areas of 100 and -100 meters high set: '100|0:0,5:5,0:10;-100|10:10,10:15,20:15'.");
            prm.declare_entry("Grid number of data points",
                              "1:1",
                              Patterns::Anything(),
                              "The number of data points of the uniform or nonuniform grid in the longitude:latitude direction.");
            // Data files
            prm.declare_entry ("Data directory",
                               "$ASPECT_SOURCE_DIR/data/geometry_model/",
                               Patterns::DirectoryName (),
                               "The name of a directory that contains the model data. ");
            prm.declare_entry ("Topography file name", "topo.dat",
                               Patterns::Anything (),
                               "The name of the file containing the topography values. In case a uniform grid is used, the file should only contain "
                               "a list of topography values for the long,lat coordinates (first longitude ascends, then latitude). A nonuniform grid "
                               "allows for variable spacing in each coordinate direction, but requires the coordinates to be specified in the file "
                               "as well. In this case, first list the longitude coordinates, then the latitude coordinates and lastly the topography values. ");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    EllipsoidalChunk<dim>::parse_parameters(ParameterHandler &prm)
    {
      // only implemented for the 3d case
      AssertThrow (dim == 3, ExcMessage ("2D has not been implemented."));
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Ellipsoidal chunk");
        {
          /**
           * Get latitude and longitudes defining region of interest from
           * the parameter file.
           */
          corners.resize(4);
          std::string NEcorner = prm.get("NE corner");
          std::string NWcorner = prm.get("NW corner");
          std::string SWcorner = prm.get("SW corner");
          std::string SEcorner = prm.get("SE corner");

          /**
           * make a list of what corners are present and check that there should be one or two corners missing,
           * otherwise throw an exception
           */
          std::vector<bool> present(4,true);
          unsigned int missing = 0;
          if (NEcorner == "")
            {
              missing++;
              present[0] = false;
            }
          if (NWcorner == "")
            {
              missing++;
              present[1] = false;
            }
          if (SWcorner == "")
            {
              missing++;
              present[2] = false;
            }
          if (SEcorner == "")
            {
              missing++;
              present[3] = false;
            }

          AssertThrow (missing != 0,
                       ExcMessage ("Only two or three of the four corners points should be provided."));
          AssertThrow (missing == 1 || missing == 2,
                       ExcMessage ("Please provide two or three corners points."));

          std::vector<double> temp;

          if (present[0])
            {
              temp = Utilities::string_to_double(Utilities::split_string_list(NEcorner,':'));
              AssertThrow (temp.size() == 2, ExcMessage ("Number of coordinates given for the NE-corner should be two (logitude:latitude)."));
              corners[0] = Point<2>(temp[0],temp[1]);
            }
          else
            corners[0]= Point<2>();

          if (present[1])
            {
              temp = Utilities::string_to_double(Utilities::split_string_list(NWcorner,':'));
              AssertThrow (temp.size() == 2, ExcMessage ("Number of coordinates given for the NW-corner should be two (logitude:latitude)."));
              corners[1] = Point<2>(temp[0],temp[1]);
            }
          else
            corners[1]= Point<2>();

          if (present[2])
            {
              temp = Utilities::string_to_double(Utilities::split_string_list(SWcorner,':'));
              AssertThrow (temp.size() == 2, ExcMessage ("Number of coordinates given for the SW-corner should be two (logitude:latitude)."));
              corners[2] = Point<2>(temp[0],temp[1]);
            }
          else
            corners[2]= Point<2>();

          if (present[3])
            {
              temp = Utilities::string_to_double(Utilities::split_string_list(SEcorner,':'));
              AssertThrow (temp.size() == 2, ExcMessage ("Number of coordinates given for the SE-corner should be two (logitude:latitude)."));
              corners[3] = Point<2>(temp[0],temp[1]);
            }
          else
            corners[3]= Point<2>();


          bottom_depth = prm.get_double("Depth");
          semi_major_axis_a = prm.get_double("Semi-major axis");
          eccentricity = prm.get_double("Eccentricity");
          semi_minor_axis_b=std::sqrt((1 - pow(eccentricity,2)) * pow(semi_major_axis_a,2));
          EW_subdiv = prm.get_double("East-West subdivisions");
          NS_subdiv = prm.get_double("North-South subdivisions");
          depth_subdiv = prm.get_double("Depth subdivisions");

          // Check whether the corners of the rectangle are really placed correctly
          if (present[0] == true && present[1] == true)
            {
              AssertThrow (corners[0][0] >= corners[1][0],
                           ExcMessage ("The longitude of the NE corner (" + boost::lexical_cast<std::string>(corners[0][0]) + ") cannot be smaller than the longitude of the NW corner (" + boost::lexical_cast<std::string>(corners[1][0]) + ")."));
            }
          if (present[0] == true && present[3] == true)
            {
              AssertThrow (corners[0][1] >= corners[3][1],
                           ExcMessage ("The latitude of the NE (" + boost::lexical_cast<std::string>(corners[0][1]) + ") corner cannot be larger than the longitude of the SE corner (" + boost::lexical_cast<std::string>(corners[3][1]) + ")."));
            }
          if (present[2] == true && present[3] == true)
            {
              AssertThrow (corners[2][0] <= corners[3][0],
                           ExcMessage ("The longitude of the SW (" + boost::lexical_cast<std::string>(corners[2][0]) + ") corner cannot be larger than the longitude of the SE corner (" + boost::lexical_cast<std::string>(corners[3][0]) + ")."));
            }
          if (present[2] == true && present[1] == true)
            {
              AssertThrow (corners[2][1] <= corners[1][1],
                           ExcMessage ("The latitude of the SW corner (" + boost::lexical_cast<std::string>(corners[2][1]) + ") cannot be smaller than the longitude of the NW corner (" + boost::lexical_cast<std::string>(corners[1][1]) + ")."));
            }
          if (missing == 2)
            {
              AssertThrow ((present[0] == true && present[2] == true) || (present[1] == true && present[3] == true),
                           ExcMessage ("Please provide two opposing corners."));

              if (present[0] == true && present[2] == true)
                AssertThrow (corners[0][0] > corners[2][0] && corners[0][1] > corners[2][1],
                             ExcMessage ("The North-East corner (" + boost::lexical_cast<std::string>(corners[0][0]) + ":"  + boost::lexical_cast<std::string>(corners[0][1]) + ") must be stricly North and East of the South-West corner (" + boost::lexical_cast<std::string>(corners[2][0]) + ":"  + boost::lexical_cast<std::string>(corners[2][1]) + ") when only two points are given."));

              if (present[1] == true && present[3] == true)
                AssertThrow (corners[1][0] < corners[3][0] && corners[1][1] > corners[3][1],
                             ExcMessage ("The North-West corner (" + boost::lexical_cast<std::string>(corners[1][0]) + ":"  + boost::lexical_cast<std::string>(corners[1][1]) + ") must be stricly North and West of the South-East corner (" + boost::lexical_cast<std::string>(corners[3][0]) + ":"  + boost::lexical_cast<std::string>(corners[3][1]) + ") when only two points are given."));
            }


          // If one or two of the corners is not provided, calculate it.
          if (missing == 1)
            {
              for (unsigned int i = 0; i <= 3; i++)
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
                  ExcMessage ("The latitude of the NE (" + boost::lexical_cast<std::string>(corners[0][1]) + ") corner cannot be larger than the longitude of the SE corner (" + boost::lexical_cast<std::string>(corners[3][1]) + "). This is an internal check, if you see this please contact the developer."));

          Assert (corners[2][0] <= corners[3][0],
                  ExcMessage ("The longitude of the SW (" + boost::lexical_cast<std::string>(corners[2][0]) + ") corner cannot be larger than the longitude of the SE corner (" + boost::lexical_cast<std::string>(corners[3][0]) + "). This is an internal check, if you see this please contact the developer."));

          Assert (corners[2][1] <= corners[1][1],
                  ExcMessage ("The latitude of the SW corner (" + boost::lexical_cast<std::string>(corners[2][1]) + ") cannot be smaller than the longitude of the NW corner (" + boost::lexical_cast<std::string>(corners[1][1]) + "). This is an internal check, if you see this please contact the developer."));

          westLongitude = corners[2][0];
          eastLongitude = corners[0][0];
          northLatitude = corners[0][1];
          southLatitude = corners[2][1];

          prm.enter_subsection("Topography");
          {
            /**
             * Determine what type of topography we have
             */
            std::string topography_type_string = prm.get("Topography type");
            TopoTypes topo_type;
            if (topography_type_string == "No topography")
              topo_type = NO_TOPOGRAPHY;
            else if (topography_type_string == "From prm exact")
              topo_type = PRM_EXACT;
            else if (topography_type_string == "From prm interpolated on uniform grid")
              topo_type = PRM_UNIFORM_GRID_INTERPOLATED;
            else if (topography_type_string == "From file uniform grid")
              topo_type = FILE_UNIFORM_GRID;
            else if (topography_type_string == "From file nonuniform grid")
              topo_type = FILE_NONUNIFORM_GRID;
            else
              {
                Assert(false,ExcMessage ("The given topography function (" + topography_type_string + ") has not been implemented."));
                topo_type = NO_TOPOGRAPHY;
              }

            manifold.topography.set_topography_type(topo_type);

            if (topo_type != NO_TOPOGRAPHY)
              {
                manifold.topography.set_corners(corners);
              }

            if (topo_type == PRM_EXACT || topo_type == PRM_UNIFORM_GRID_INTERPOLATED)
              {
                /**
                 * we need to fill the point lists and topography values. They
                 * are stored in the Topography subsection in the Topography parameter.
                 */
                std::vector<double> topography_values;
                std::vector<std::vector<std::vector<double> > > point_lists;
                std::string temptopo = prm.get("Topography parameters");
                std::vector<std::string> temp_topographies = Utilities::split_string_list(temptopo,';');
                unsigned int temp_topographies_size = temp_topographies.size();

                topography_values.resize(temp_topographies_size,0);
                point_lists.resize(temp_topographies_size);
                for (unsigned int i_topo = 0; i_topo < temp_topographies_size; i_topo++)
                  {
                    std::vector<std::string> temp_topography = Utilities::split_string_list(temp_topographies[i_topo],'|');
                    Assert(temp_topography.size() == 2,ExcMessage ("The given line '" + temp_topographies[i_topo] + "' is not correct. "
                                                                   "It consists of " + boost::lexical_cast<std::string>(temp_topography.size())
                                                                   + " parts separated by |, but it should only contain 2 parts: "
                                                                   "the height and the coordinate, separated by a |.'"));

                    topography_values[i_topo] = Utilities::string_to_double(temp_topography[0]);

                    std::vector<std::string> temp_coordinates = Utilities::split_string_list(temp_topography[1],',');
                    unsigned int temp_coordinate_size = temp_coordinates.size();
                    point_lists[i_topo].resize(temp_coordinate_size,std::vector<double>(2,0));

                    for (unsigned int i_coord = 0; i_coord < temp_coordinate_size; i_coord++)
                      {
                        point_lists[i_topo][i_coord] = Utilities::string_to_double(Utilities::split_string_list(temp_coordinates[i_coord],':'));
                        Assert(temp_topography.size() == 2,ExcMessage ("The given coordinate '" + temp_coordinates[i_coord] + "' is not correct. "
                                                                       "It consists of " + boost::lexical_cast<std::string>(temp_topography.size())
                                                                       + " parts separated by :, but it should only contain 2 parts: "
                                                                       "the longitude and latitude of the the coordinate, separated by a ':'."));
                      }

                    manifold.topography.set_topography_values (topography_values);
                    manifold.topography.set_point_lists (point_lists);
                  }
              }

            if (topo_type == PRM_UNIFORM_GRID_INTERPOLATED || topo_type == FILE_UNIFORM_GRID || topo_type == FILE_NONUNIFORM_GRID)
              {
                std::vector<double> grid_number_data_points = Utilities::string_to_double(Utilities::split_string_list(prm.get("Grid number of data points"),':'));
                AssertThrow(grid_number_data_points.size() == 2, ExcMessage("The number of grid points needs to be specified for the longitude and latitude directions. "));
                manifold.topography.set_grid_number_data_points(grid_number_data_points);
              }

            if (topo_type == FILE_UNIFORM_GRID || topo_type == FILE_NONUNIFORM_GRID)
              {
                std::string data_directory  = prm.get ("Data directory");
                {
                  const std::string      subst_text = "$ASPECT_SOURCE_DIR";
                  std::string::size_type position;
                  while (position = data_directory.find (subst_text),  position!=std::string::npos)
                    data_directory.replace (data_directory.begin()+position,
                                            data_directory.begin()+position+subst_text.size(),
                                            subst_text);
                }
                const std::string topo_file_name = prm.get("Topography file name");

                manifold.topography.set_topography_file(data_directory+topo_file_name);
              }

          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      /**
       * Construct manifold object Pointer to an object that describes the geometry.
       */
      manifold.set_manifold_parameters(semi_major_axis_a,
                                       eccentricity,
                                       semi_minor_axis_b,
                                       bottom_depth,
                                       corners);
    }

    template <int dim>
    double
    EllipsoidalChunk<dim>::depth(const Point<dim> &position) const
    {
      return std::max(std::min(-manifold.pull_back(position)[2], maximal_depth()), 0.0);
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
      const Point<dim> long_lat_depth = manifold.pull_back(position);
      return semi_major_axis_a / (std::sqrt(1 - eccentricity * eccentricity * std::sin(long_lat_depth[1]) * std::sin(long_lat_depth[1])));
    }

    template <int dim>
    double
    EllipsoidalChunk<dim>::get_eccentricity() const
    {
      return eccentricity;
    }

    template <int dim>
    const std::vector<Point<2> > &
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
      Point<dim> p =
        (manifold.push_forward(Point<3>(southLatitude,
                                        eastLongitude,
                                        -bottom_depth))
         + manifold.push_forward(Point<3>(northLatitude, eastLongitude, 0)))
        / 2;
      p /= p.norm();
      p *= get_radius(p) - depth;
      return p;
    }
  }
}

template <int dim>
typename aspect::GeometryModel::EllipsoidalChunk<dim>::EllipsoidalChunkGeometry
aspect::GeometryModel::EllipsoidalChunk<dim>::get_manifold() const
{
  return manifold;
}

// explicit instantiations
namespace aspect
{
  namespace GeometryModel
  {
    ASPECT_REGISTER_GEOMETRY_MODEL(EllipsoidalChunk,
                                   "ellipsoidal chunk",
                                   "A 3D chunk geometry that accounts for Earth's ellipticity (default assuming the "
                                   "WGS84 ellipsoid definition) which can be defined in non-coordinate directions. "
                                   "In the description of the ellipsoidal chunk, two of the ellipsoidal axes have the "
                                   "same length so that there is only a semi-major axis and a semi-minor axis. "
                                   "The user has two options for creating an ellipsoidal chunk geometry: 1) by defining "
                                   "two opposing points (SW and NE or NW and SE) a coordinate parallel ellipsoidal "
                                   "chunk geometry will be created. 2) by defining three points a non-coordinate "
                                   "parallel ellipsoidal chunk will be created. The points are defined in the input "
                                   "file by longitude:latitude. It is also possible to define additional subdivisions of the "
                                   "mesh in each direction. Faces of the model are defined as 0, west; 1,east; 2, south; 3, "
                                   "north; 4, inner; 5, outer. "
                                   "Additionally, one can add topography to the ellipsoidal chunk, either directly from "
                                   "the input file, or by reading in a data file. ")
  }
}
