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
    double vof_from_d<2> (Tensor<1, 2, double> normal,
                          double d)
    {
      const int dim = 2;
      double norm1, max, mpos, dtest;

      //Get 1-Norm
      norm1 = 0.0;
      max = 0.0;
      for (unsigned int i = 0; i < dim; ++i)
        {
          double term = numbers::NumberTraits<double>::abs (normal[i]);
          norm1 += term;
          max = (max < term) ? term : max;
        }

      //Obtain volume
      if (d <= -0.5*norm1)
        {
          return 0.0;
        }
      if (d >= 0.5*norm1)
        {
          return 1.0;
        }
      dtest = d / norm1;
      mpos = 1.0 - max/norm1;
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
    double d_from_vof<2> (Tensor<1, 2, double> normal,
                          double vol)
    {
      const int dim = 2;
      double norm1, max, mpos;

      //Get 1-Norm
      norm1 = 0.0;
      max = 0.0;
      for (unsigned int i = 0; i < dim; ++i)
        {
          double term = numbers::NumberTraits<double>::abs (normal[i]);
          norm1 += term;
          max = (max < term) ? term : max;
        }

      if (norm1 == 0.0)
        {
          norm1 = 1.0;
          mpos = 0.0;
        }
      else
        {
          mpos = 1.0 - max/norm1;
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
    double vof_from_d<3> (Tensor<1, 3, double> normal,
                          double d)
    {
      // 3D vol calculation not yet implemented
      // Likely to require additional calculations
      return 0.0;
    }

    template<>
    double d_from_vof<3> (Tensor<1, 3, double> normal,
                          double vol)
    {
      // 3D interface location calculation not yet implemented
      // Almost certain to require iterative method
      return 0.0;
    }

    template<int dim>
    double calc_vof_flux_edge (Tensor<1, dim, double> dir,
                               Tensor<1, dim, double> vflux,
                               Tensor<1, dim, double> normal,
                               double d)
    {
      // If flux outward, return 0
      if (dir*vflux<=0.0)
        {
          return 0.0;
        }
      Tensor<1, dim, double> i_normal;
      double i_d;

      i_d = d+0.5*(normal*vflux-dir*normal);
      for (unsigned int i = 0; i<dim; ++i)
        {
          if (dir[i] == 0.0)
            {
              i_normal[i] = normal[i];
            }
          else
            {
              i_normal[i] = normal*vflux;
            }
        }
      return vof_from_d (i_normal, i_d);
    }

    template
    double calc_vof_flux_edge<2>(Tensor<1, 2, double> dir,
                                 Tensor<1, 2, double> vflux,
                                 Tensor<1, 2, double> normal,
                                 double d);
    template
    double calc_vof_flux_edge<3>(Tensor<1, 3, double> dir,
                                 Tensor<1, 3, double> vflux,
                                 Tensor<1, 3, double> normal,
                                 double d);

  }
}
