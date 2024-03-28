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

#include "world_builder/objects/natural_coordinate.h"

namespace WorldBuilder
{
  namespace Objects
  {
    NaturalCoordinate::NaturalCoordinate(const std::array<double,3> &position,
                                         const CoordinateSystems::Interface &coordinate_system_)
    {
      coordinate_system = coordinate_system_.natural_coordinate_system();
      coordinates = coordinate_system_.cartesian_to_natural_coordinates(position);
    }


    // todo, should be possible to make this without the interface, since the Point knows the coord system.
    NaturalCoordinate::NaturalCoordinate(const Point<3> &position,
                                         const CoordinateSystems::Interface &coordinate_system_)
    {
      coordinate_system = coordinate_system_.natural_coordinate_system();
      coordinates = coordinate_system_.cartesian_to_natural_coordinates(position.get_array());
    }


    const std::array<double,3> &NaturalCoordinate::get_coordinates() const
    {
      return coordinates;
    }


    std::array<double,2> NaturalCoordinate::get_surface_coordinates() const
    {
      std::array<double,2> coordinate;

      switch (coordinate_system)
        {
          case CoordinateSystem::cartesian:
            coordinate[0] = coordinates[0];
            coordinate[1] = coordinates[1];
            break;

          case CoordinateSystem::spherical:
            coordinate[0] = coordinates[1];
            coordinate[1] = coordinates[2];
            break;

          default:
            WBAssertThrow (false, "Coordinate system not implemented.");
        }

      return coordinate;
    }


    Point<2> NaturalCoordinate::get_surface_point() const
    {
      Point<2> coordinate(0,0,coordinate_system);

      switch (coordinate_system)
        {
          case CoordinateSystem::cartesian:
            coordinate[0] = coordinates[0];
            coordinate[1] = coordinates[1];
            break;

          case CoordinateSystem::spherical:
            coordinate[0] = coordinates[1];
            coordinate[1] = coordinates[2];
            break;

          default:
            WBAssertThrow (false, "Coordinate system not implemented.");
        }

      return coordinate;
    }


    CoordinateSystem
    NaturalCoordinate::get_coordinate_system() const
    {
      return coordinate_system;
    }


    double &NaturalCoordinate::get_ref_depth_coordinate()
    {
      switch (coordinate_system)
        {
          case CoordinateSystem::cartesian:
            return coordinates[2];

          case CoordinateSystem::spherical:
            return coordinates[0];

          default:
            WBAssertThrow (false, "Coordinate system not implemented.");
        }

      return coordinates[2];
    }
  } // namespace Objects
} // namespace WorldBuilder