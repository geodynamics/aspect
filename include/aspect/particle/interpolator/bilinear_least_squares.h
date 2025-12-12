/*
 Copyright (C) 2017 - 2024 by the authors of the ASPECT code.

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

#ifndef _aspect_particle_interpolator_bilinear_least_squares_h
#define _aspect_particle_interpolator_bilinear_least_squares_h

#include <aspect/particle/interpolator/linear_least_squares.h>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      /**
       * The LinearLeastSquares interpolator used to use a bilinear least
       * squares algorithm and was called BilinearLeastSquares. Since the
       * algorithm was changed the name had to change as well. The new algorithm
       * is marginally less accurate, but faster and less prone to
       * oscillations. 
       * 
       * @deprecated This name is deprecated and will be removed.
       */
      DEAL_II_DEPRECATED_WITH_COMMENT("The class <BilinearLeastSquares> is now named "
        "<LinearLeastSquares> and will be removed in the future.")
      using BilinearLeastSquares = LinearLeastSquares;
    }
  }
}

#endif
