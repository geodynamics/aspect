/*
  Copyright (C) 2018-2024 by the authors of the World Builder code.

  This file is part of the World Builder.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "world_builder/coordinate_systems/spherical.h"


#include "world_builder/types/double.h"
#include "world_builder/types/object.h"
#include "world_builder/utilities.h"

namespace WorldBuilder
{
  namespace CoordinateSystems
  {
    Spherical::Spherical(WorldBuilder::World *world_)
    {
      this->world = world_;
    }

    Spherical::~Spherical()
      = default;

    void
    Spherical::declare_entries(Parameters &prm, const std::string & /*unused*/)
    {

      // Add depth method to the required parameters.
      prm.declare_entry("", Types::Object({"depth method"}),
                        "A spherical coordinate system. The coordinates are (radius, longitude, latitude). "
                        "The radius is set in this plugin, the longitude extends at least from -360 to 360 degrees, "
                        "and the latitude extends from -90 to 90. It is required to choose a depth method. Please "
                        "see the manual for more information.");


      prm.declare_entry("depth method",
                        Types::String("",std::vector<std::string>({"starting point", "begin segment","begin at end segment", "continuous"})),
                        R"(Which depth method to use in the spherical case. The available options are 'starting point', )"
                        R"('begin segment' and 'begin at end segment'. See the manual section on coordinate systems for )"
                        R"(more info.)");

      prm.declare_entry("radius",
                        Types::Double(6371000.),
                        R"(The radius of the sphere.)");


    }

    void
    Spherical::parse_entries(Parameters &prm)
    {
      prm.enter_subsection("coordinate system");
      {
        const std::string string_depth_method = prm.get<std::string>("depth method");
        if (string_depth_method == "starting point")
          used_depth_method = DepthMethod::angle_at_starting_point_with_surface;
        else if (string_depth_method == "begin segment")
          used_depth_method = DepthMethod::angle_at_begin_segment_with_surface;
        else if (string_depth_method == "begin at end segment")
          used_depth_method = DepthMethod::angle_at_begin_segment_applied_to_end_segment_with_surface;
        //else if (string_depth_method == "continuous")
        //used_depth_method = DepthMethod::continuous_angle_with_surface;
        else
          WBAssertThrow(true,"Option " << string_depth_method << " is not a valid depth method for spherical "
                        "coordinates. The available options are 'starting point', 'begin segment' and 'begin at end segment'. "
                        "The option 'continuous' is not yet available.");

        radius_sphere = prm.get<double>("radius");
      }
      prm.leave_subsection();
    }


    CoordinateSystem
    Spherical::natural_coordinate_system() const
    {
      return CoordinateSystem::spherical;
    }


    DepthMethod
    Spherical::depth_method() const
    {
      return used_depth_method;
    }


    std::array<double,3>
    Spherical::cartesian_to_natural_coordinates(const std::array<double,3> &position) const
    {
      return Utilities::cartesian_to_spherical_coordinates(Point<3>(position,cartesian));
    }


    std::array<double,3>
    Spherical::natural_to_cartesian_coordinates(const std::array<double,3> &position) const
    {
      return Utilities::spherical_to_cartesian_coordinates(position).get_array();
    }

    double
    Spherical::distance_between_points_at_same_depth(const Point<3> &point_1, const Point<3> &point_2) const
    {
      WBAssert(point_1.get_coordinate_system() == spherical,
               "Can not convert non spherical points through the spherical coordinate system.");
      WBAssert(point_2.get_coordinate_system() == spherical,
               "Can not convert non spherical points through the spherical coordinate system.");
      const double radius = point_1[0];
      WBAssert((radius - point_2[0]) < std::numeric_limits<double>::epsilon() * std::max(1.0,radius), "The radius of point 1 is not the same as the radius of point 2.");

      // based on https://math.stackexchange.com/questions/1304169/distance-between-two-points-on-a-sphere
      const Point<3> point_1_cart = Utilities::spherical_to_cartesian_coordinates(point_1.get_array());
      const Point<3> point_2_cart = Utilities::spherical_to_cartesian_coordinates(point_2.get_array());

      return radius * std::acos(std::min(1.,std::max(0.,point_1_cart*point_2_cart/(radius*radius))));
    }


    double
    Spherical::max_model_depth() const
    {
      return radius_sphere;
    }

    /**
     * Register plugin
     */
    WB_REGISTER_COORDINATE_SYSTEM(Spherical, spherical)
  } // namespace CoordinateSystems
} // namespace WorldBuilder

