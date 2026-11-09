// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2025-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include <aspect/material_model/thermal_conductivity/constant.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {
      template <int dim>
      void
      Constant<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &/*in*/,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        for (auto &thermal_conductivity: out.thermal_conductivities)
          thermal_conductivity = k;
      }



      template <int dim>
      void
      Constant<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Thermal conductivity", "4.7",
                           Patterns::Double (0.),
                           "The value of the thermal conductivity $k$. "
                           "Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");
      }



      template <int dim>
      void
      Constant<dim>::parse_parameters (ParameterHandler &prm)
      {
        k = prm.get_double ("Thermal conductivity");
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {
#define INSTANTIATE(dim) \
  template class Constant<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}
