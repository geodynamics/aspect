/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


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
