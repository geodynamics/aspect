/*
  Copyright (C) 2016 by the authors of the ASPECT code.

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

#include <aspect/vofinterface/VOFUtils.h>

namespace aspect
{
  namespace InterfaceTracker
  {
    using namespace dealii;

    template<>
    double vol_from_d<2> (Tensor<1, 2, double> normal,
                          double d)
    {
      const int dim = 2;
      double norm1, mpos, dtest;

      //Get 1-Norm
      norm1 = 0.0;
      for (unsigned int i = 0; i < dim; ++i)
        {
          norm1 += numbers::NumberTraits<double>::abs (normal[i]);
        }

      //Find min component
      mpos = 0.6;
      for (unsigned int i = 0; i < dim; ++i)
        {
          double mcand = numbers::NumberTraits<double>::abs (normal[i])
                         / norm1;
          mpos = (mcand < mpos) ? mcand : mpos;
        }

      //Obtain volume
      dtest = d / norm1;
      if (dtest <= -0.5)
        {
          return 0.0;
        }
      if (dtest >= 0.5)
        {
          return 1.0;
        }
      if (dtest < mpos - 0.5)
        {
          return (dtest + 0.5) * (dtest + 0.5) / (2.0*mpos * (1.0 - mpos));
        }
      if (dtest > 0.5 - mpos)
        {
          return 1.0 - (dtest - 0.5) * (dtest - 0.5) / (2.0*mpos * (1.0 - mpos));
        }
      return 0.5 + dtest / (1.0 - mpos);
    }

    template<>
    double d_from_vol<2> (Tensor<1, 2, double> normal,
                          double vol)
    {
      const int dim = 2;
      double norm1, mpos;

      //Get 1-Norm
      norm1 = 0.0;
      for (unsigned int i = 0; i < dim; ++i)
        {
          norm1 += numbers::NumberTraits<double>::abs (normal[i]);
        }

      //Find min component
      mpos = 0.6;
      for (unsigned int i = 0; i < dim; ++i)
        {
          double mcand = numbers::NumberTraits<double>::abs (normal[i])
                         / norm1;
          mpos = (mcand < mpos) ? mcand : mpos;
        }

      //Obtain const
      if (vol <= 0.0)
        {
          return -0.5 * norm1;
        }
      if (vol >= 1.0)
        {
          return 0.5 * norm1;
        }
      if (vol < 0.5 * mpos / (1 - mpos))
        {
          return norm1 * (-0.5 + sqrt (2.0*vol * mpos * (1 - mpos)));
        }
      if (vol > 1.0 - 0.5 * mpos / (1 - mpos))
        {
          return norm1 * (0.5 - sqrt (2.0*(1.0 - vol) * mpos * (1 - mpos)));
        }
      return norm1 * (1 - mpos) * (vol - 0.5);

    }

    template<>
    double vol_from_d<3> (Tensor<1, 3, double> normal,
                          double d)
    {
      // 3D vol calculation not yet implemented
      // Likely to require additional calculations
      return 0.0;
    }

    template<>
    double d_from_vol<3> (Tensor<1, 3, double> normal,
                          double vol)
    {
      // 3D interface location calculation not yet implemented
      // Almost certain to require iterative method
      return 0.0;
    }

  }
}
