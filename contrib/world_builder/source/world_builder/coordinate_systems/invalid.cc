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

#include "world_builder/coordinate_systems/invalid.h"
#include "world_builder/types/object.h"


namespace WorldBuilder
{
  namespace CoordinateSystems
  {
    Invalid::Invalid(WorldBuilder::World *world_)
    {
      this->world = world_;
    }

    Invalid::~Invalid()
      = default;

    void
    Invalid::parse_entries(Parameters & /*prm*/)
    {}


    CoordinateSystem
    Invalid::natural_coordinate_system() const
    {
      return CoordinateSystem::invalid;
    }


    DepthMethod
    Invalid::depth_method() const
    {
      return DepthMethod::none;
    }


    std::array<double,3>
    Invalid::cartesian_to_natural_coordinates(const std::array<double,3> & /*position*/) const
    {
      return {{NaN::DQNAN,NaN::DQNAN,NaN::DQNAN}};
    }


    std::array<double,3>
    Invalid::natural_to_cartesian_coordinates(const std::array<double,3> & /*position*/) const
    {
      return {{NaN::DQNAN,NaN::DQNAN,NaN::DQNAN}};
    }


    double
    Invalid::distance_between_points_at_same_depth(const Point<3> & /*point_1*/, const Point<3> & /*point_2*/) const
    {
      return NaN::DQNAN;
    }


    double
    Invalid::max_model_depth() const
    {
      return NaN::DQNAN;
    }

    /**
     * Not registering plugin because this is only used for testing.
     */
    //WB_REGISTER_COORDINATE_SYSTEM(Invalid, invalid)
  } // namespace CoordinateSystems
} // namespace WorldBuilder

