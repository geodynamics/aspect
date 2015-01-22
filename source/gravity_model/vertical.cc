/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include <aspect/gravity_model/vertical.h>

#include <deal.II/base/tensor.h>

namespace aspect
{
  namespace GravityModel
  {
    template <int dim>
    Tensor<1,dim>
    Vertical<dim>::gravity_vector (const Point<dim> &) const
    {
      Tensor<1,dim> g;
      g[dim-1] = -gravity_magnitude;
      return g;
    }
    template <int dim>
    void
    Vertical<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Vertical");
        {
          prm.declare_entry ("Magnitude", "1",
                             Patterns::Double (0),
                             "Value of the gravity vector in $m/s^2$ directed "
                             "along negative y (2D) or z (3D) axis.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    Vertical<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Vertical");
        {
          gravity_magnitude = prm.get_double ("Magnitude");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace GravityModel
  {
    ASPECT_REGISTER_GRAVITY_MODEL(Vertical,
                                  "vertical",
                                  "A gravity model in which the gravity direction is vertically downward "
                                  "and at a constant magnitude by default equal to one.")
  }
}
