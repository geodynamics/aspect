/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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


#include <aspect/gravity_model/radial_constant.h>

#include <aspect/geometry_model/interface.h>
#include <aspect/coordinate_systems.h>

#include <deal.II/base/tensor.h>

namespace aspect
{
  namespace GravityModel
  {
    template <int dim>
    Tensor<1,dim>
    RadialConstant<dim>::gravity_vector (const Point<dim> &p) const
    {
      if (p.norm() == 0.0)
        return Tensor<1,dim>();

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
                             Patterns::Double (),
                             "Magnitude of the gravity vector in $m/s^2$. For positive values "
                             "the direction is radially inward towards the center of the earth.");
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

      AssertThrow (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical ||
                   this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::ellipsoidal,
                   ExcMessage ("Gravity model 'radial constant' should not be used with geometry models that "
                               "do not have either a spherical or ellipsoidal natural coordinate system."));
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
                                  "A gravity model in which the gravity has a constant magnitude "
                                  "and the direction is radial (pointing inward if the value "
                                  "is positive). The magnitude is read from the parameter "
                                  "file in subsection 'Radial constant'.")
  }
}
