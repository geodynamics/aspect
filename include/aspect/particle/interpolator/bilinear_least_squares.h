// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2017-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

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
      template <int dim>
      using BilinearLeastSquares
      DEAL_II_DEPRECATED_WITH_COMMENT("The class <BilinearLeastSquares> is now named "
                                      "<LinearLeastSquares> and will be removed in the future.")
        = LinearLeastSquares<dim>;
    }
  }
}

#endif
