/*
  Copyright (C) 2019 by the authors of the ASPECT code.

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


#include <aspect/material_model/rheology/constant_viscosity.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      double
      ConstantViscosity::compute_viscosity () const
      {
        return eta;
      }



      void
      ConstantViscosity::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Viscosity", "1e21",
                           Patterns::Double (0),
                           "The value of the viscosity $\\eta$. Units: $kg/m/s$ or $Pa s$.");
      }



      void
      ConstantViscosity::parse_parameters (ParameterHandler &prm)
      {
        eta = prm.get_double ("Viscosity");
      }
    }
  }
}
