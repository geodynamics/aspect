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

        // COMPUTE NATURAL LOGARITHM OF PRESSURE AND TEMPERATURE
        // ==========================================================================
          P_log = std::log(in.pressures[i]);
          T_log = std::log(in.temperatures[i]);
        // ==========================================================================

        // DEFINE COEFFICIENTS FOR LATTICE THERMAL CONDUCTIVITY OF DIFFERENT MINERALS
        // ==========================================================================
        // Coefficients for Dry Olivine
          OlivineDry_LatTC_a0 =  -4.124100000;
          OlivineDry_LatTC_b1 =   2.146900000;
          OlivineDry_LatTC_ymin = 1.280933845;
          OlivineDry_LatTC_ymax = 2.607263820;
        // ==========================================================================
        // Coefficients for Orthopyroxene (Enstatite)
          OpxEnstati_LatTC_a0 =   -3.004700000;
          OpxEnstati_LatTC_b1 =    2.600000000;
          OpxEnstati_LatTC_ymin =  1.760865151; 
          OpxEnstati_LatTC_ymax =  2.096937429;
        // ==========================================================================
        // Coefficients for Garnet (Pyrope)
          GrtPyropes_LatTC_a0 =   -4.363700000;
          GrtPyropes_LatTC_b1 =    2.036800000;
          GrtPyropes_LatTC_ymin =  1.481604541; 
          GrtPyropes_LatTC_ymax =  2.443131606;
        // ==========================================================================

        // COMPUTE LATTICE THERMAL CONDUCTIVITIES IN REAL AND SIMPLEX SPACE
        // ==========================================================================
        // Latice Thermal Conductivity of Dry Olivine
          OlivineDry_LatTC_zSimpl = OlivineDry_LatTC_a0 + OlivineDry_LatTC_b1*P_log;
          OlivineDry_LatTC_ySimpl = std::exp(OlivineDry_LatTC_zSimpl);
          OlivineDry_LatTC_yPrime = OlivineDry_LatTC_ySimpl/(1+OlivineDry_LatTC_ySimpl);
          OlivineDry_LatTC_yTCLat = OlivineDry_LatTC_ymin+(OlivineDry_LatTC_ymax-OlivineDry_LatTC_ymin)*OlivineDry_LatTC_yPrime;
          OlivineDry_Lattice_TCon = std::exp(OlivineDry_LatTC_yTCLat);
        // ==========================================================================
        // Latice Thermal Conductivity of Orthopyroxene (Enstatite)
          OpxEnstati_LatTC_zSimpl = OpxEnstati_LatTC_a0 + OpxEnstati_LatTC_b1*P_log;
          OpxEnstati_LatTC_ySimpl = std::exp(OpxEnstati_LatTC_zSimpl);
          OpxEnstati_LatTC_yPrime = OpxEnstati_LatTC_ySimpl/(1+OpxEnstati_LatTC_ySimpl);
          OpxEnstati_LatTC_yTCLat = OpxEnstati_LatTC_ymin+(OpxEnstati_LatTC_ymax-OpxEnstati_LatTC_ymin)*OpxEnstati_LatTC_yPrime;
          OpxEnstati_Lattice_TCon = std::exp(OpxEnstati_LatTC_yTCLat);
          // ==========================================================================
        // Latice Thermal Conductivity of Garnet (Pyrope)
          GrtPyropes_LatTC_zSimpl = GrtPyropes_LatTC_a0 + GrtPyropes_LatTC_b1*P_log;
          GrtPyropes_LatTC_ySimpl = std::exp(GrtPyropes_LatTC_zSimpl);
          GrtPyropes_LatTC_yPrime = GrtPyropes_LatTC_ySimpl/(1+GrtPyropes_LatTC_ySimpl);
          GrtPyropes_LatTC_yTCLat = GrtPyropes_LatTC_ymin+(GrtPyropes_LatTC_ymax-GrtPyropes_LatTC_ymin)*GrtPyropes_LatTC_yPrime;
          GrtPyropes_Lattice_TCon = std::exp(GrtPyropes_LatTC_yTCLat);
          // ==========================================================================

          out.thermal_conductivities[i] = OlivineDry_Lattice_TCon;

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