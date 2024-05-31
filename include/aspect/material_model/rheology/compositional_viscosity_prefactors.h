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
    using namespace dealii;

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

          /**
           * Compute the viscosity.
           */
          double
          compute_viscosity (const MaterialModel::MaterialModelInputs<dim> &in,
                             const double base_viscosity,
                             const unsigned int composition_index,
                             const unsigned int q) const;

        private:
          /**
           * The viscosity prefactors or terms used to calculate the viscosity
           * prefactors, which are read in from the input file by the
           * parse_parameters() function. Users can choose between different schemes.
           * none: no prefactor multiplication
           * hk04_olivine_hydration: calculate the water fugacity from
           * Hirth & Kohlstaedt 2004 10.1029/138GM06. using the compositional field
           * 'bound_fluid'.
           * The prefactor for a given compositional field is multiplied with a
           * base_viscosity value provided by the material model, which
           * is then returned to the material model.
           */
          enum ViscosityPrefactorScheme
          {
            none,
            hk04_olivine_hydration
          } viscosity_prefactor_scheme;

          // Initialize variables for the water fugacity calculation, from HK04
          std::vector<double> water_fugacity_exponents;
          // From Hirth & Kohlstaedt 2004, equation 6
          const double activation_energy_H2O = 40e3;
          const double activation_volume_H2O = 10e-6;

          // We calculate the Molar mass of olivine using the molar mass of fayalite (203.79) and the
          // molar mass of forsterite (140.693), and a mass fraction of 90% forsterite.
          const double M_olivine = 147.0027;
          const double M_H2O = 18.01528; // Molar mass of H2O
      };
    }
  }
}
#endif
