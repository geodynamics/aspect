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
        
        // DEFINE VOLUME FRACTIONS OF MINERALS
        // ==========================================================================
        const double OlivineDry_fraction   =  0.6;
        const double OpxEnstati_fraction   =  0.2;
        const double GrtPyropes_fraction   =  0.2;
        // ==========================================================================

        // DEFINE COEFFICIENTS FOR LATTICE THERMAL CONDUCTIVITY OF DIFFERENT MINERALS
        // ==========================================================================
        // Coefficients for Dry Olivine
        const double OlivineDry_LatTC_a0 =   -4.124100000;
        const double OlivineDry_LatTC_b1 =    2.146900000;
        const double OlivineDry_LatTC_ymin =  1.280933845;
        const double OlivineDry_LatTC_ymax =  2.607263820;
        const double OlivineDry_TDepnd_Exp = -0.5;      
        // ==========================================================================
        // Coefficients for Orthopyroxene (Enstatite)
        const double OpxEnstati_LatTC_a0 =   -3.004700000;
        const double OpxEnstati_LatTC_b1 =    2.600000000;
        const double OpxEnstati_LatTC_ymin =  1.760865151; 
        const double OpxEnstati_LatTC_ymax =  2.096937429;
        const double OpxEnstati_TDepnd_Exp = -0.5;
        // ==========================================================================
        // Coefficients for Garnet (Pyrope)
        const double GrtPyropes_LatTC_a0 =   -4.363700000;
        const double GrtPyropes_LatTC_b1 =    2.036800000;
        const double GrtPyropes_LatTC_ymin =  1.481604541; 
        const double GrtPyropes_LatTC_ymax =  2.443131606;
        const double GrtPyropes_TDepnd_Exp = -0.4314;
        // ==========================================================================

        // DEFINE COEFFICIENTS FOR RADIATIVE THERMAL CONDUCTIVITY OF DIFFERENT MINERALS
        // ==========================================================================
        // Coefficients for Dry Olivine
        const double OlivineDry_RadTC_a0 =   -10.00900000;
        const double OlivineDry_RadTC_b1 =    1.883900000;
        const double OlivineDry_RadTC_ymin = -23.02585093;
        const double OlivineDry_RadTC_ymax =  1.289885976;
        // ==========================================================================
        // Coefficients for Orthopyroxene (Enstatite)
        const double OpxEnstati_RadTC_a0 =   -13.532000000;
        const double OpxEnstati_RadTC_b1 =    2.4004000000;
        const double OpxEnstati_RadTC_ymin = -23.025850930; 
        const double OpxEnstati_RadTC_ymax =  1.4456685920;
        // ==========================================================================
        // Coefficients for Garnet (Pyrope)
        const double GrtPyropes_RadTC_a0 =   -11.782000000;
        const double GrtPyropes_RadTC_b1 =    2.0718000000;
        const double GrtPyropes_RadTC_ymin = -23.025850930; 
        const double GrtPyropes_RadTC_ymax =  1.4479836950;
        // ==========================================================================

        // DEFINE ROOM TEMPERATURE [K] AND MODEL TEMPERATURE
        // ==========================================================================
        const double T_room = 298.15; 
        // ==========================================================================

        for (unsigned int i = 0; i < n_points; ++i) 
        {

        // COMPUTE NATURAL LOGARITHM OF PRESSURE AND TEMPERATURE
        // ==========================================================================
          double P_log = std::log(in.pressures[i]);
          double T_log = std::log(in.temperatures[i]);
        // ==========================================================================

        // TAKE MODEL TEMPERATURE
        // ==========================================================================
          double T_mod = in.temperatures[i];
        // ==========================================================================

        // COMPUTE THE LATTICE THERMAL CONDUCTIVITY BOUNDARIES IN REAL (+,-) AND SIMPLEX (0->1) SPACE
        // ==========================================================================
          // Dry Olivine
          double OlivineDry_LatTC_zSimpl = OlivineDry_LatTC_a0 + OlivineDry_LatTC_b1*P_log;
          double OlivineDry_LatTC_ySimpl = std::exp(OlivineDry_LatTC_zSimpl);
          double OlivineDry_LatTC_yPrime = OlivineDry_LatTC_ySimpl/(1+OlivineDry_LatTC_ySimpl);
          double OlivineDry_LatTC_yTCLat = OlivineDry_LatTC_ymin+(OlivineDry_LatTC_ymax-OlivineDry_LatTC_ymin)*OlivineDry_LatTC_yPrime;
        // ==========================================================================
          // Orthopyroxene (Enstatite)
          double OpxEnstati_LatTC_zSimpl = OpxEnstati_LatTC_a0 + OpxEnstati_LatTC_b1*P_log;
          double OpxEnstati_LatTC_ySimpl = std::exp(OpxEnstati_LatTC_zSimpl);
          double OpxEnstati_LatTC_yPrime = OpxEnstati_LatTC_ySimpl/(1+OpxEnstati_LatTC_ySimpl);
          double OpxEnstati_LatTC_yTCLat = OpxEnstati_LatTC_ymin+(OpxEnstati_LatTC_ymax-OpxEnstati_LatTC_ymin)*OpxEnstati_LatTC_yPrime;
        // ==========================================================================
          // Garnet (Pyrope)
          double GrtPyropes_LatTC_zSimpl = GrtPyropes_LatTC_a0 + GrtPyropes_LatTC_b1*P_log;
          double GrtPyropes_LatTC_ySimpl = std::exp(GrtPyropes_LatTC_zSimpl);
          double GrtPyropes_LatTC_yPrime = GrtPyropes_LatTC_ySimpl/(1+GrtPyropes_LatTC_ySimpl);
          double GrtPyropes_LatTC_yTCLat = GrtPyropes_LatTC_ymin+(GrtPyropes_LatTC_ymax-GrtPyropes_LatTC_ymin)*GrtPyropes_LatTC_yPrime;
        // ==========================================================================

        // COMPUTE P-DEPENDENT LATTICE THERMAL CONDUCTIVITY OF MINERALS
        // ==========================================================================
          // Dry Olivine
          double OlivineDry_PDep_LatTCon = std::exp(OlivineDry_LatTC_yTCLat);
          // Orthopyroxene (Enstatite)
          double OpxEnstati_PDep_LatTCon = std::exp(OpxEnstati_LatTC_yTCLat);
          // Garnet (Pyrope)
          double GrtPyropes_PDep_LatTCon = std::exp(GrtPyropes_LatTC_yTCLat);
        // ==========================================================================

        // COMPUTE T-DEPENDENT LATTICE THERMAL CONDUCTIVITY OF MINERALS
        // ==========================================================================
          // Dry Olivine
          double OlivineDry_TDep_LatTCon = std::pow(OlivineDry_PDep_LatTCon*(T_mod/T_room),OlivineDry_TDepnd_Exp);
          // Orthopyroxene (Enstatite)
          double OpxEnstati_TDep_LatTCon = std::pow(OpxEnstati_PDep_LatTCon*(T_mod/T_room),OpxEnstati_TDepnd_Exp);
          // Garnet (Pyrope)
          double GrtPyropes_TDep_LatTCon = std::pow(GrtPyropes_PDep_LatTCon*(T_mod/T_room),GrtPyropes_TDepnd_Exp);
        // ==========================================================================

        // COMPUTE P,T-DEPENDENT LATTICE THERMAL CONDUCTIVITY OF MINERALS
        // ==========================================================================
          // Dry Olivine
          double OlivineDry_PTDep_LatTCo = OlivineDry_TDep_LatTCon;
          // Orthopyroxene (Enstatite)
          double OpxEnstati_PTDep_LatTCo = OpxEnstati_TDep_LatTCon;
          // Garnet (Pyrope)
          double GrtPyropes_PTDep_LatTCo = GrtPyropes_TDep_LatTCon;
        // ==========================================================================


          
  

            out.thermal_conductivities[i] = OlivineDry_PTDep_LatTCo;
          }

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