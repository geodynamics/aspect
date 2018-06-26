/*
  Copyright (C) 2018 by the authors of the World Builder code.

  This file is part of the World Builder.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef _aspect_world_builder_coordinate_sytems_h
#define _aspect_world_builder_coordinate_sytems_h

  namespace WorldBuilder
  {
        /**
         * This enum lists available coordinate systems that can be used for
         * the function variables. Allowed values are 'cartesian',
         * 'spherical', and 'depth'. 'spherical' coordinates follow: r, phi
         * (2D) or r, phi, theta (3D); where r is radius, phi is longitude,
         * and theta is the polar angle (colatitude).
         */
        enum CoordinateSystem
        {
          cartesian,
          spherical,
          invalid
        };
      }

#endif

