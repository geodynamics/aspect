// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_geometry_model_initial_topography_model_zero_topography_h
#define _aspect_geometry_model_initial_topography_model_zero_topography_h

#include <aspect/geometry_model/initial_topography_model/interface.h>

namespace aspect
{
  namespace InitialTopographyModel
  {
    /**
     * A class that implements zero initial topography.
     *
     * @ingroup InitialTopographyModels
     */
    template <int dim>
    class ZeroTopography : public Interface<dim>
    {
      public:
        /**
         * Return the value of the initial topography as a function of position.
         *
         * For the current class, this function obviously simply returns a zero
         * value.
         */
        double
        value (const Point<dim-1> &p) const override;

        /**
         * Return the maximum value of the elevation.
         */
        double max_topography () const override;
    };
  }
}


#endif
