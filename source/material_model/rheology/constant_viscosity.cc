// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2019-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/material_model/rheology/constant_viscosity.h>

#include <deal.II/base/signaling_nan.h>


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
