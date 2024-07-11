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

#include "world_builder/coordinate_systems/cartesian.h"
#include "world_builder/types/object.h"


namespace WorldBuilder
{
  namespace CoordinateSystems
  {
    Cartesian::Cartesian(WorldBuilder::World *world_)
    {
      this->world = world_;
    }

    Cartesian::~Cartesian()
      = default;

    void
    Cartesian::declare_entries(Parameters &prm, const std::string & /*unused*/)
    {
      prm.declare_entry("", Types::Object(),
                        "A Cartesian coordinate system. Coordinates are (x,y,z) and extend infinitely in all directions.");
    }

    void
    Cartesian::parse_entries(Parameters & /*prm*/)
    {}


    CoordinateSystem
    Cartesian::natural_coordinate_system() const
    {
      return CoordinateSystem::cartesian;
    }


    DepthMethod
    Cartesian::depth_method() const
    {
      return DepthMethod::none;
    }


    std::array<double,3>
    Cartesian::cartesian_to_natural_coordinates(const std::array<double,3> &position) const
    {
      return position;
    }


    std::array<double,3>
    Cartesian::natural_to_cartesian_coordinates(const std::array<double,3> &position) const
    {
      return position;
    }


    double
    Cartesian::distance_between_points_at_same_depth(const Point<3> &point_1, const Point<3> &point_2) const
    {
      WBAssert(point_1.get_coordinate_system() == cartesian,
               "Can not convert non-Cartesian points through the Cartesian coordinate system.");
      WBAssert(point_2.get_coordinate_system() == cartesian,
               "Can not convert non-Cartesian points through the Cartesian coordinate system.");
      // Todo: check that points are at the same depth.
      const Point<3> difference = point_1-point_2;
      const Point<2> point_at_depth(difference[0],difference[1], cartesian);

      return point_at_depth.norm();
    }


    double
    Cartesian::max_model_depth() const
    {
      return std::numeric_limits<double>::infinity();
    }

    /**
     * Register plugin
     */
    WB_REGISTER_COORDINATE_SYSTEM(Cartesian, cartesian)
  } // namespace CoordinateSystems
} // namespace WorldBuilder

