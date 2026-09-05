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
           * Return @p base_viscosity multiplied by the viscosity factor for
           * composition @p composition_index at evaluation point @p q. For
           * the Peng-Robinson scheme, the viscosity is calculated as
           * $\eta=\eta_{\mathrm{base}}f_{\mathrm{H_2O}}^{-r/n}$. The fugacity
           * $f_{\mathrm{H_2O}}$ depends only on temperature and pressure and
           * is therefore independent of composition, although $r/n$ may be
           * configured separately for each compositional field. The
           * @p modified_flow_laws argument selects the exponent for diffusion
           * or dislocation creep.
           */
          double
          compute_viscosity (const MaterialModel::MaterialModelInputs<dim> &in,
                             const double base_viscosity,
                             const unsigned int composition_index,
                             const unsigned int q,
                             const ModifiedFlowLaws &modified_flow_laws) const;

          /**
           * Compute the fugacity of the configured pure fluid from the
           * Peng-Robinson equation of state at @p temperature in K and
           * absolute @p pressure in Pa.
           * The calculation is independent of composition: at a given
           * temperature and pressure, every composition has the same fugacity.
           * The returned fugacity is in Pa. An exception is raised if the pressure
           * exceeds the equation of state limit of 2.5 GPa.
           * Fluid-saturated conditions are assumed.
           */
          double
          compute_fugacity (const double temperature, const double pressure) const;

          /**
           * Return whether the selected viscosity prefactor scheme computes
           * fugacity with the Peng-Robinson equation of state.
           */
          bool
          uses_peng_robinson_fugacity () const;

          /**
           * Create the named Peng-Robinson fugacity output when the
           * corresponding viscosity-prefactor scheme is selected.
           */
          void
          create_fugacity_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const;

          /**
           * Fill the named Peng-Robinson fugacity output at evaluation point
           * @p point_index, if the output object was requested.
           */
          void
          fill_fugacity_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
                                 const unsigned int point_index,
                                 MaterialModel::MaterialModelOutputs<dim> &out) const;

        private:
          /**
           * The viscosity prefactors or terms used to calculate the viscosity
           * prefactors, which are read in from the input file by the
           * parse_parameters() function. Users can choose between different schemes.
           * none: no viscosity change
           * hk04_olivine_hydration: calculate the viscosity change due to hydrogen
           * incorporation into olivine using Hirth & Kohlstedt 2004 10.1029/138GM06.
           * This method requires a composition called 'bound_fluid' which tracks the wt%
           * water in the solid, which is used to compute an atomic ratio of H/Si ppm
           * assuming 90 mol% forsterite and 10 mol% fayalite, and finally calculates
           * a water fugacity.
           * peng_robinson76_fugacity: estimate a viscosity modification
           * using the fugacity of a pure fluid at a given temperature and
           * pressure determined using the
           * Peng-Robinson equation of state introduced by Peng & Robinson
           * (1976, 10.1021/i160057a011). The viscosity modification
           * assumes that the material is saturated in the pure fluid. The viscosity is
           * calculated as
           * $\eta=\eta_{\mathrm{base}}f_{\mathrm{H_2O}}^{-r/n}$. Fugacity is
           * composition-independent because it depends only on temperature
           * and pressure, although $r/n$ may be different for
           * each compositional field. This formulation approximates the experimentally
           * observed weakening of minerals by water (e.g., Karato et al., 1986,
           * 10.1029/JB091iB08p08151), with the water-fugacity term in the
           * creep flow law following Mei & Kohlstedt (2000,
           * 10.1029/2000JB900180).
           * interface weakening: calculate the viscosity change due to the presence of
           * a sub-grid scale weak layer present at the interface of two other compositions.
           */
          enum ViscosityPrefactorScheme
          {
            none,
            hk04_olivine_hydration,
            peng_robinson76_fugacity,
            interface_weakening,
          };
          /**
           *  This variable is read from the parameter file through a parameter called 'Viscosity prefactor scheme'.
           */
          ViscosityPrefactorScheme viscosity_prefactor_scheme;

          // Initialize variables for the water fugacity calculation, from HK04
          /**
           *  This variable is read from the parameter file through a parameter called 'Water fugacity exponents for diffusion creep'.
           */
          std::vector<double> diffusion_water_fugacity_exponents;
          /**
           *  This variable is read from the parameter file through a parameter called 'Water fugacity exponents for dislocation creep'.
           */
          std::vector<double> dislocation_water_fugacity_exponents;
          /**
           *  This variable is read from the parameter file through a parameter called 'Minimum mass fraction bound water content for fugacity'.
           */
          std::vector<double> minimum_mass_fraction_water_for_dry_creep;

          // Variables for the interface weakening scheme
          std::vector<std::string> weakening_field_names;
          std::vector<double> interface_weakening_factors;
          double interface_weakening_threshold;

          // From Hirth & Kohlstedt 2004, equation 6
          const double A_H2O = 2.6e-5; // 1/Pa
          const double activation_energy_H2O = 40e3; // J/mol/K
          const double activation_volume_H2O = 10e-6; // m^3/mol

          // We calculate the molar mass of olivine using the molar mass of fayalite (0.20379 kg/mol)
          // and the molar mass of forsterite (0.140693 kg/mol), and a mole fraction of 90% forsterite
          // in olivine.
          const double molar_mass_olivine = 0.1470027; // kg/mol
          const double molar_mass_H2O = 0.01801528; // kg/mol

          /**
           * Critical temperature of the fluid in K for the Peng-Robinson
           * equation of state. This value is fluid-specific and user defined.
           */
          double critical_temperature;

          /**
           * Critical pressure of the fluid in Pa for the Peng-Robinson
           * equation of state. This value is fluid-specific and user defined.
           */
          double critical_pressure;

          /**
           * Peng-Robinson substance constant calculated from the acentric
           * factor using equation (18) of Peng & Robinson (1976,
           * 10.1021/i160057a011).
           */
          double kappa;

          /**
           * Acentric factor of the fluid. Peng & Robinson describe kappa as
           * a constant characteristic of each substance and correlate it
           * against the acentric factor in equation (18) of Peng & Robinson
           * (1976, 10.1021/i160057a011).
           */
          double acentric_factor;

          /**
           * Maximum pressure in Pa used in the Peng-Robinson fugacity calculation.
           * Pressures above this value are set to the cutoff.
           */
          double pressure_cutoff;

      };
    }
  }
}
#endif
