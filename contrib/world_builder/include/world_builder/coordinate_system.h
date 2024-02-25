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

#ifndef WORLD_BUILDER_COORDINATE_SYSTEM_H
#define WORLD_BUILDER_COORDINATE_SYSTEM_H

namespace WorldBuilder
{
  /**
   * This enum lists available coordinate systems that can be used for the
   * function variables. Allowed values are 'cartesian', 'spherical'.
   * 'spherical' coordinates follow: r, phi (2D) or r, phi, theta (3D); where
   * r is radius, phi is longitude, and theta is the polar angle (colatitude).
   */
  enum CoordinateSystem
  {
    cartesian,
    spherical,
    invalid
  };

  /**
   * There are three ways to look at what should happen when the user gives an
   * angle to go down in a sphere. The first one is the easiest, which is just
   * to tread the system as a Cartesian system. This means taking the angle
   * with the surface tangential at the starting point at the surface
   * everywhere. The advantage is that it is very intuitive and easy to
   * predict what geometry should be produced. This option is called
   * 'angle_at_starting_point_with_surface'. The other end member is that
   * going down with an angle in a sphere should mean that it always goes down
   * with a constant angle relative to the surface tangential at that
   * longitude and latitude. This means a the geometry would a a logarithmic
   * spiral. This option is called 'continuous_angle_with_surface'. Because it
   * is currently to time consuming to implement this option with all the
   * required features, there is an option which approximates this behavior,
   * by using the option 'angle_at_starting_point_with_surface'. This option
   * uses the surface tangential at the beginning of every segment. This means
   * that in the limit of segments going to length zero, this option produces
   * exactly the same geometry as 'continuous_angle_with_surface'.
   */
  enum DepthMethod
  {
    angle_at_starting_point_with_surface,
    angle_at_begin_segment_with_surface,
    angle_at_begin_segment_applied_to_end_segment_with_surface,
    continuous_angle_with_surface,
    none
  };
} // namespace WorldBuilder

#endif

