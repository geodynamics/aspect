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

#ifndef WORLD_BUILDER_OBJECTS_NATURAL_COORDINATE_H
#define WORLD_BUILDER_OBJECTS_NATURAL_COORDINATE_H

#include <array>
#include "world_builder/point.h"
#include "world_builder/coordinate_systems/interface.h"

namespace WorldBuilder
{
  namespace Objects
  {
    /*
    * A class that represents a point in a chosen coordinate system.
    */
    class NaturalCoordinate
    {
      public:
        /**
         * Constructor based on providing the geometry model as a pointer
         */
        NaturalCoordinate(const std::array<double,3> &position,
                          const ::WorldBuilder::CoordinateSystems::Interface &coordinate_system);

        /**
         * Constructor based on providing the geometry model as a pointer
         */
        NaturalCoordinate(const Point<3> &position,
                          const ::WorldBuilder::CoordinateSystems::Interface &coordinate_system);

        /**
         * Returns the coordinates in the given coordinate system, which may
         * not be Cartesian.
         */
        const std::array<double,3> &get_coordinates() const;

        /**
         * The coordinate that represents the 'surface' directions in the
         * chosen coordinate system.
         */
        std::array<double,2> get_surface_coordinates() const;

        /**
         * The coordinate that represents the 'surface' directions in the
         * chosen coordinate system.
         */
        Point<2> get_surface_point() const;

        /**
         * The coordinate that represents the 'depth' direction in the chosen
         * coordinate system.
         */
        double get_depth_coordinate() const;

        /**
         * Return a reference to the coordinate that represents the 'depth' direction
         * in the chosen coordinate system.
         */
        double &get_ref_depth_coordinate();

        /**
         * get the coordinate system type of this coordinate.
         */
        CoordinateSystem get_coordinate_system() const;

      private:
        /**
         * An enum which stores the coordinate system of this natural
         * point
         */
        CoordinateSystem coordinate_system;

        /**
         * An array which stores the coordinates in the coordinates system
         */
        std::array<double,3> coordinates;
    };


    inline
    double
    NaturalCoordinate::get_depth_coordinate() const
    {
      switch (coordinate_system)
        {
          case CoordinateSystem::cartesian:
          {
            return coordinates[2];
            break;
          }

          case CoordinateSystem::spherical:
          {
            return coordinates[0];
            break;
          }

          default:
            WBAssertThrow (false, "Coordinate system not implemented.");
        }

      return 0;
    }
  }
}

#endif