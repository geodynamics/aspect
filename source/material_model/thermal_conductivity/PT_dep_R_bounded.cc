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

#include <aspect/material_model/thermal_conductivity/PT_dep_R_bounded.h>

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

        // Coefficients for Dry Olivine
          OlivineDry_a0 =  -4.124100000;
          OlivineDry_b1 =   2.146900000;
          OlivineDry_ymin = 1.280933845;
          OlivineDry_ymax = 2.607263820;

        // Compute natural logarithm of pressure
          P_log = std::log(in.pressures[i]);

        // Compute thermal conductivity
          OlivineDry_zSimpl = OlivineDry_a0 + OlivineDry_b1*P_log;
          OlivineDry_ySimpl = std::exp(OlivineDry_zSimpl);
          OlivineDry_yPrime = OlivineDry_ySimpl/(1+OlivineDry_ySimpl);
          OlivineDry_yTCLat = OlivineDry_ymin+(OlivineDry_ymax-OlivineDry_ymin)*OlivineDry_yPrime;
        
          out.thermal_conductivities[i] = std::exp(OlivineDry_yTCLat);
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