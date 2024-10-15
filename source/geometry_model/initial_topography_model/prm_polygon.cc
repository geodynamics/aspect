/*
  Copyright (C) 2016 - 2024 by the authors of the ASPECT code.

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


#include <aspect/geometry_model/initial_topography_model/prm_polygon.h>
#include <aspect/utilities.h>
#include <boost/lexical_cast.hpp>

namespace aspect
{
  namespace InitialTopographyModel
  {
    template <int dim>
    double
    PrmPolygon<dim>::
    value (const Point<dim-1> &p) const
    {
      const Point<2> p1 = (dim == 2 ? Point<2>(p[0],0) : Point<2>(p[0],p[1]));

      /**
       * We go through the loop in the reverse order, because we
       * want the last entry in the list to prescribe the value.
       */
      for (unsigned int i = point_lists.size(); i > 0; i--)
        {
          if (aspect::Utilities::polygon_contains_point<dim>(point_lists[i-1],p1))
            {
              return topography_values[i-1];
            }
        }

      // Point is not in any of the polygons. Return zero.
      return 0;
    }



    template <int dim>
    double
    PrmPolygon<dim>::
    max_topography () const
    {
      return maximum_topography;
    }



    template <int dim>
    void
    PrmPolygon<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Initial topography model");
        {
          prm.enter_subsection("Prm polygon");
          {
            prm.declare_entry("Topography parameters",
                              "",
                              Patterns::Anything(),
                              "Set the topography height and the polygon which should be set to that height. "
                              "The format is : \"The topography height \textgreater The point list describing "
                              "a polygon \\& The next topography height \textgreater the next point list "
                              "describing a polygon.\" The format for the point list describing the polygon is "
                              "\"x1,y1;x2,y2\". For example for two triangular areas of 100 and -100 meters high "
                              "set: '100 \textgreater 0,0;5,5;0,10 \\& -100 \textgreater 10,10;10,15;20,15'. "
                              "Units of the height are always in meters. The units of the coordinates are "
                              "dependent on the geometry model. In the box model they are in meters, in the "
                              "chunks they are in degrees, etc. Please refer to the manual of the individual "
                              "geometry model to so see how the topography is implemented.");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    PrmPolygon<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Initial topography model");
        {
          prm.enter_subsection("Prm polygon");
          {
            /**
             * we need to fill the point lists and topography values. They
             * are stored in the Topography subsection in the Topography parameter.
             */
            maximum_topography = std::numeric_limits<double>::lowest();
            const std::string temptopo = prm.get("Topography parameters");
            const std::vector<std::string> temp_topographies = Utilities::split_string_list(temptopo,'&');
            const unsigned int temp_topographies_size = temp_topographies.size();

            topography_values.resize(temp_topographies_size,0);
            point_lists.resize(temp_topographies_size);
            for (unsigned int i_topo = 0; i_topo < temp_topographies_size; ++i_topo)
              {
                const std::vector<std::string> temp_topography = Utilities::split_string_list(temp_topographies[i_topo],'>');
                Assert(temp_topography.size() == 2,ExcMessage ("The given line '" + temp_topographies[i_topo]
                                                               + "' is not correct. It consists of "
                                                               + boost::lexical_cast<std::string>(temp_topography.size())
                                                               + " parts separated by >, but it should only contain "
                                                               "two parts: the height and the list of point coordinates,"
                                                               " separated by a >.'"));

                topography_values[i_topo] = Utilities::string_to_double(temp_topography[0]);
                maximum_topography = std::max(topography_values[i_topo],maximum_topography);

                const std::vector<std::string> temp_coordinates = Utilities::split_string_list(temp_topography[1],';');
                const unsigned int temp_coordinate_size = temp_coordinates.size();
                point_lists[i_topo].resize(temp_coordinate_size);
                for (unsigned int i_coord = 0; i_coord < temp_coordinate_size; ++i_coord)
                  {
                    const std::vector<double> temp_point = Utilities::string_to_double(Utilities::split_string_list(temp_coordinates[i_coord],','));
                    Assert(temp_point.size() == 2,ExcMessage ("The given coordinate '" + temp_coordinates[i_coord] + "' is not correct. "
                                                              "It consists of " + boost::lexical_cast<std::string>(temp_topography.size())
                                                              + " parts separated by :, but it should only contain 2 parts: "
                                                              "the two coordinates of the polygon points, separated by a ','."));


                    point_lists[i_topo][i_coord] = Point<2>(temp_point[0], temp_point[1]);
                  }

              }
          }
          prm.leave_subsection();
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
  namespace InitialTopographyModel
  {
    ASPECT_REGISTER_INITIAL_TOPOGRAPHY_MODEL(PrmPolygon,
                                             "prm polygon",
                                             "An initial topography model that defines the initial topography "
                                             "as constant inside each of a set of polygonal parts of the "
                                             "surface. The polygons, and their associated surface elevation, "
                                             "are defined in the `Geometry model/Initial topography/Prm polygon' "
                                             "section.")
  }
}
