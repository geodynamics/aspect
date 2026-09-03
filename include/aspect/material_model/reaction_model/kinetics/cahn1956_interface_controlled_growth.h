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

#ifndef _aspect_material_model_reaction_model_kinetics_cahn1956_interface_controlled_growth_h
#define _aspect_material_model_reaction_model_kinetics_cahn1956_interface_controlled_growth_h

#include <aspect/material_model/reaction_model/kinetics/cahn1956_interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      /**
       * A kinetic model that computes the reaction rate of a phase transformation governed by `interface-controlled growth` kinetics. This
       * describes phase transformations where a single product phase nucleates along grain boundaries of a single reactant phase, and then
       * grows via grain-boundary migration at the expense of the reactant phase. The kinetic model is derived from the site-saturated
       * framework of Cahn (1956) and with an interface-controlled growth rate after Fisher (1978; their Equation (11) ).
       *
       * Examples of this type of reaction include the olivine -> wadsleyite transformation (Hosoya et al., 2005), aragonite -> calcite,
       * wadsleyite -> ringwoodite and coesite -> quartz. In principle, the same functional form is applicable to any transformation for
       * which diffusionless grain-boundary migration is rate-limiting.
       *
       * The net forward reaction rate dX_B/dt accounts for both the forward (A -> B) and reverse (B -> A) transformations, where X_A = 1 -
       * X_B:
       *
       *                  ┌  Z * T * exp(-(Ha+(P*Va))/(R*T)) * (1-exp(dG/(R*T))) * (1-X_B) if dG < 0 (dX_B/dt is positive, A -> B)
       *      dX_B/dt =  ─┤
       *                  └  Z * T * exp(-(Ha+(P*Va))/(R*T)) * (1-exp(dG/(R*T))) * X_B if dG > 0 (dX_B/dt is negative, B -> A)
       *
       * where X_B is the cumulative reaction progress (volume fraction of B; V_B / (V_A + V_B)), Z is a kinetic prefactor, Ha and Va are
       * the activation enthalpy and volume of the reaction, and dG is the Gibbs energy change of the forward reaction (dG < 0 favors
       * product B; dG > 0 favors reactant A).
       *
       * This plugin can store parameters for multiple reactions.
       *
       * The formulation above comes from Equations (1) and (2) of Hosoya et al. (2005), with modified notation to match the
       * EutectoidDecomposition kinetic model (after Kerswell et al., 2026).
       *
       * Cahn, J. W. (1956). https://doi.org/10.1016/0001-6160(56)90041-4
       * Fisher, G. W. (1978). https://doi.org/10.1016/0016-7037(78)90292-2
       * Hosoya et al. (2002). https://doi.org/10.1029/2005GL023398
       * Kerswell et al. (2026). https://doi.org/10.1029/2026JB033781
       *
       * @ingroup ReactionModel
       */
      template <int dim>
      class InterfaceControlledGrowth : public Cahn1956Interface<dim>
      {
        public:
          /**
           * Compute the net forward reaction rate dX_B/dt (units: 1/s, or 1/yr if the ``Use years instead of seconds'' global parameter is
           * set),
           * for the reaction at local index @p reaction_index, given the local temperature (K), pressure (Pa), Gibbs energy change of the forward
           * reaction (A -> B; J/mol), and the cumulative forward reaction progress field from the Material Model.
           *
           * The sign of dX_B/dt is determined by delta_gibbs_energy:
           * - If dG < 0 (favoring product B), the returned rate is positive (dX_B/dt > 0) and limited by X_A = 1 - X_B
           * - If dG > 0 (favoring reactant A), the returned rate is negative (dX_B/dt < 0) and limited by X_B
           */
          double
          net_forward_reaction_rate(const double       temperature,
                                    const double       pressure,
                                    const double       delta_forward_gibbs_energy,
                                    const double       cumulative_forward_reaction_progress,
                                    const unsigned int reaction_index) const override;

          /**
           * Declare the parameters from the parameter file.
           */
          static void
          declare_parameters(ParameterHandler &prm);

          /**
           * Read the parameters from the parameter file.
           */
          void
          parse_parameters(ParameterHandler &prm, const unsigned int n_reactions) override;

        private:
          /**
           * Kinetic prefactor Z, one per reaction using this model. Units: 1/s/K.
           */
          std::vector<double> kinetic_prefactors;

          /**
           * Activation enthalpy Ha of the reaction, one per reaction using this model. Units: J/mol.
           */
          std::vector<double> activation_enthalpies;

          /**
           * Activation volume Va of the reaction, one per reaction using this model. Units: m^3/mol.
           */
          std::vector<double> activation_volumes;
      };
    }
  }
}

#endif
