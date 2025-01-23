/*
  Copyright (C) 2024 - by the authors of the ASPECT code.

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

#include <aspect/material_model/thermal_conductivity/constant.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {
      template <int dim>
      void
      Constant<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        const unsigned int n_points = in.n_evaluation_points();
        for (unsigned int i = 0; i < n_points; ++i)
          out.thermal_conductivities[i] = k;
      }



      template <int dim>
      void
      Constant<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Material model");
        {
          prm.enter_subsection("Thermal conductivity");
          {
            prm.enter_subsection("Constant");
            {
              prm.declare_entry ("Thermal conductivity", "4.7",
                                 Patterns::Double (0.),
                                 "The value of the thermal conductivity $k$. "
                                 "Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }



      template <int dim>
      void
      Constant<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Material model");
        {
          prm.enter_subsection("Thermal conductivity");
          {
            prm.enter_subsection("Constant");
            {
              k = prm.get_double ("Thermal conductivity");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
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
