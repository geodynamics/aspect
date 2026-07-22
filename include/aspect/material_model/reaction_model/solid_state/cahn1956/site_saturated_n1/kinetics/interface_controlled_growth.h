/*
  Copyright (C) 2011 - 2026 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_reaction_model_solid_state_cahn1956_site_saturated_n1_kinetics_interface_controlled_growth_h
#define _aspect_material_model_reaction_model_solid_state_cahn1956_site_saturated_n1_kinetics_interface_controlled_growth_h

#include <aspect/material_model/reaction_model/solid_state/cahn1956/site_saturated_n1/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      namespace SolidState
      {
        namespace Cahn1956
        {
          namespace SiteSaturatedN1
          {
            /**
             * A reaction model that computes the rate of a polymorphic (solid-solid) phase transformation under interface-controlled growth kinetics,
             * following the formulation of Hosoya et al. (2005) for the olivine -> wadsleyite transformation. The same functional form is applicable to
             * any two-phase transformation for which nucleation and interface migration, rather than long-range diffusion, are rate-limiting.
             *
             * The reaction rate for the forward transformation A -> B is:
             *     dX/dt = Z * T * exp( -(Ha + P * Va) / (R * T) ) * (1 - exp(dG / (R * T))) * (1 - X_B)
             * where X_A = (1 - X_B) is the mass (or volume) fraction of the reactant phase A, Z is a kinetic prefactor, Ha and Va are the activation enthalpy
             * and volume of the interface reaction, and dG is the Gibbs free energy change of the transformation (dG < 0 favors the product phase B).
             *
             * This class is agnostic to how dG is computed (from a thermodynamic data table, analytic Clausius-Clapeyron approximation, etc.) and
             * agnostic to how the resulting rate is used (compositional field reaction term, latent heat, viscosity weakening, etc.). Caller supplies
             * dG and X_A. For the reverse reaction (when dG > 0), remember to pass the mass fraction of the product phase instead. You may need to
             * enforce phase-fraction constraints (e.g., mass balancing).
             *
             * @ingroup ReactionModel
             */
            template <int dim>
            class InterfaceControlledGrowth : public Interface<dim>
            {
              public:
                /**
                 * Declare the parameters this function takes through input files.
                 */
                static void declare_parameters(ParameterHandler &prm);

                /**
                 * Read the parameters from the parameter file.
                 */
                void parse_parameters(ParameterHandler &prm) override;

                /**
                 * Compute the reaction rate dX/dt (units: 1/s, or 1/yr if the ``Use years instead of seconds'' global parameter is set) for the forward
                 * transformation A -> B, given the local temperature (K), pressure (Pa), the Gibbs free energy change delta_gibbs_energy (J/mol) of the
                 * transformation, and the mass or volume fraction reactant_phase_fraction of the phase being consumed.
                 *
                 * The sign convention follows delta_gibbs_energy: if it is negative (phase B is thermodynamically favored), the returned rate is
                 * positive and proportional to reactant_phase_fraction (interpreted as X_A, the fraction of phase A). If it is positive (phase A is
                 * favored), the returned rate is negative and reactant_phase_fraction should instead be supplied as X_B, the fraction of phase B.
                 */
                double reaction_rate(const double temperature, const double pressure, const double delta_gibbs_energy, const double reactant_phase_fraction) const override;

                /**
                 * Return the Arrhenius factor exp( -(Ha + P * Va) / (R * T) ) on its own.
                 */
                double arrhenius_factor(const double temperature, const double pressure) const override;

                /**
                 * Return the (signed) thermodynamic driving force (1 - exp(-|dG| / (R * T))) on its own.
                 */
                double thermodynamic_factor(const double temperature, const double delta_gibbs_energy) const override;

              private:
                /**
                 * Kinetic prefactor Z. Units: 1/s/K.
                 */
                double kinetic_factor;

                /**
                 * Activation enthalpy Ha of the interface reaction. Units: J/mol.
                 */
                double activation_enthalpy;

                /**
                 * Activation volume Va of the interface reaction. Units: m^3/mol.
                 */
                double activation_volume;

                /**
                 * Universal gas constant. Units: J/mol/K.
                 */
                static constexpr double gas_constant = 8.314;

                /**
                 * Evaluate exp(x), clamping the exponent to +-700 to avoid floating point overflow/underflow for pathological inputs.
                 */
                static double clamped_exp(const double x);
            };
          }
        }
      }
    }
  }
}

#endif
