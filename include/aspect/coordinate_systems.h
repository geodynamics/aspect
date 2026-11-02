// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2018-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_coordinate_systems_h
#define _aspect_coordinate_systems_h

#include <deal.II/base/utilities.h>


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
