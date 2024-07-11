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

#ifndef WORLD_BUILDER_OBJECTS_DISTANCE_FROM_SURFACE_H
#define WORLD_BUILDER_OBJECTS_DISTANCE_FROM_SURFACE_H

namespace WorldBuilder
{
  namespace Objects
  {
    class PlaneDistances
    {
      public:

        /**
         * Constructor to create an PlaneDistances object from the distance
         * from and along the surface of plane
         */
        PlaneDistances(const double distance_from_surface_, const double distance_along_surface_)
          :
          distance_from_surface(distance_from_surface_),
          distance_along_surface(distance_along_surface_)
        {}

        /**
         * Function to get the distance from the surface
        */
        double get_distance_from_surface() const;

        /**
         * Function to get the distance along the surface
        */
        double get_distance_along_surface() const;

      private:
        /**
        * the distance from the surface of the plane
        */
        double distance_from_surface;
        /**
        * the distance along the surface of the plane
        */
        double distance_along_surface;
    };
  }
}

#endif