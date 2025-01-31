/*
  Copyright (C) 2024 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_rheology_compositional_viscosity_prefactors_h
#define _aspect_material_model_rheology_compositional_viscosity_prefactors_h

#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      /**
       * A class that handles multiplication of viscosity for a given compositional
       * field. The multiplication factors for each composition (viscosity
       * prefactors) are also declared, parsed, and in some cases calculated in this class.
       */
      template <int dim>
      class CompositionalViscosityPrefactors : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          CompositionalViscosityPrefactors();

          /**
           * Declare the parameters this function takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm);

          // The flow laws that can be
          // currently modified.
          enum ModifiedFlowLaws
          {
            diffusion,
            dislocation
          } modified_flow_laws;

          /**
           * Compute the viscosity.
           */
          double
          compute_viscosity (const MaterialModel::MaterialModelInputs<dim> &in,
                             const double base_viscosity,
                             const unsigned int composition_index,
                             const unsigned int q,
                             const ModifiedFlowLaws &modified_flow_laws) const;

        private:
          /**
           * The viscosity prefactors or terms used to calculate the viscosity
           * prefactors, which are read in from the input file by the
           * parse_parameters() function. Users can choose between different schemes.
           * none: no viscosity change
           * hk04_olivine_hydration: calculate the viscosity change due to hydrogen
           * incorporation into olivine using Hirth & Kohlstaedt 2004 10.1029/138GM06.
           * This method requires a composition called 'bound_fluid' which tracks the wt%
           * water in the solid, which is used to compute an atomic ratio of H/Si ppm
           * assuming 90 mol% forsterite and 10 mol% fayalite, and finally calculates
           * a water fugacity.
           * The prefactor for a given compositional field is multiplied with a
           * base_viscosity value provided by the material model, which is then returned
           * to the material model.
           */
          enum ViscosityPrefactorScheme
          {
            none,
            hk04_olivine_hydration,
          } viscosity_prefactor_scheme;

          // Initialize variables for the water fugacity calculation, from HK04
          std::vector<double> diffusion_water_fugacity_exponents;
          std::vector<double> dislocation_water_fugacity_exponents;
          std::vector<double> minimum_mass_fraction_water_for_dry_creep;

          // From Hirth & Kohlstaedt 2004, equation 6
          const double A_H2O = 2.6e-5; // 1/Pa
          const double activation_energy_H2O = 40e3; // J/mol/K
          const double activation_volume_H2O = 10e-6; // m^3/mol

          // We calculate the molar mass of olivine using the molar mass of fayalite (0.20379 kg/mol)
          // and the molar mass of forsterite (0.140693 kg/mol), and a mole fraction of 90% forsterite
          // in olivine.
          const double molar_mass_olivine = 0.1470027; // kg/mol
          const double molar_mass_H2O = 0.01801528; // kg/mol
      };
    }
  }
}
#endif
