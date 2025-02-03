/*
  Copyright (C) 2018 - 2023 by the authors of the ASPECT code.

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

#ifndef _aspect_coordinate_systems_h
#define _aspect_coordinate_systems_h

namespace aspect
{
  namespace Utilities
  {
    /**
     * Because many places in ASPECT assume that all functions in the namespace
     * <code>dealii::Utilities</code> are available without qualification as
     * <code>Utilities::function</code>, just as all the function in the
     * namespace <code>aspect::Utilities</code>, we make sure all these functions
     * are available inside <code>aspect::Utilities</code>. This is maybe not
     * the cleanest solution, but it is most compatible with a lot of existing
     * code, and also allows to migrate ASPECT functions into deal.II when
     * useful without introducing incompatibilities.
     *
     * We need to do this in every header that introduces something into the
     * namespace <code>aspect::Utilities</code>, because it needs to happen
     * no matter which header files of ASPECT are included.
     */
    using namespace dealii::Utilities;

    namespace Coordinates
    {
      /**
       * This enum lists available coordinate systems that can be used for
       * the function variables. Allowed values are 'cartesian',
       * 'spherical', and 'depth'. 'spherical' coordinates follow: r, phi
       * (2D) or r, phi, theta (3D); where r is radius, phi is longitude,
       * and theta is the polar angle (colatitude). The 'depth' is a
       * one-dimensional coordinate system in which only the distance
       * below the 'top' surface (depth) as defined by each geometry model,
       * is used.
       */
      enum CoordinateSystem
      {
        depth,
        cartesian,
        spherical,
        ellipsoidal,
        invalid
      };
    }
  }
}

#endif
