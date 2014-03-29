/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
/*  $Id$  */


#include <aspect/gravity_model/radial.h>

#include <deal.II/base/tensor.h>

namespace aspect
{
  namespace GravityModel
  {
// ------------------------------ RadialConstant -------------------
    template <int dim>
    Tensor<1,dim>
    RadialConstant<dim>::gravity_vector (const Point<dim> &p) const
    {
      const double r = p.norm();
      return -magnitude * p/r;
    }



    template <int dim>
    void
    RadialConstant<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Radial constant");
        {
          prm.declare_entry ("Magnitude", "9.81",
                             Patterns::Double (0),
                             "Magnitude of the gravity vector in $m/s^2$. The direction is "
                             "always radially inward towards the center of the earth.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    RadialConstant<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Radial constant");
        {
          magnitude = prm.get_double ("Magnitude");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }

// ------------------------------ RadialEarthLike -------------------

    template <int dim>
    Tensor<1,dim>
    RadialEarthLike<dim>::gravity_vector (const Point<dim> &p) const
    {
      const double r = p.norm();
      return -(1.245e-6 * r + 7.714e13/r/r) * p / r;
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace GravityModel
  {
    ASPECT_REGISTER_GRAVITY_MODEL(RadialConstant,
                                  "radial constant",
                                  "A gravity model in which the gravity direction is radially inward "
                                  "and at constant magnitude. The magnitude is read from the parameter "
                                  "file in subsection 'Radial constant'.")

    ASPECT_REGISTER_GRAVITY_MODEL(RadialEarthLike,
                                  "radial earth-like",
                                  "A gravity model in which the gravity direction is radially inward "
                                  "and with a magnitude that matches that of the earth at "
                                  "the core-mantle boundary as well as at the surface and "
                                  "in between is physically correct under the assumption "
                                  "of a constant density.")
  }
}
