/*
  Copyright (C) 2019 - 2020 by the authors of the ASPECT code.

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

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      ConstantViscosity::ConstantViscosity ()
        :
        viscosity(numbers::signaling_nan<double>())
      {}



      double
      ConstantViscosity::compute_viscosity () const
      {
        return viscosity;
      }



      void
      ConstantViscosity::declare_parameters (ParameterHandler &prm,
                                             const double default_viscosity)
      {
        prm.declare_entry ("Viscosity", std::to_string(default_viscosity),
                           Patterns::Double (0.),
                           "The value of the viscosity $\\eta$. Units: \\si{\\pascal\\second}.");
      }



      void
      ConstantViscosity::parse_parameters (ParameterHandler &prm)
      {
        viscosity = prm.get_double ("Viscosity");
      }
    }
  }
}
