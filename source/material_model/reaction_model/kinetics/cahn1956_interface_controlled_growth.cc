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

#include <aspect/global.h>

#include <aspect/material_model/reaction_model/kinetics/cahn1956_interface_controlled_growth.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/patterns.h>

#include <cmath>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      template <int dim>
      double
      InterfaceControlledGrowth<dim>::net_forward_reaction_rate(const double       temperature,
                                                                const double       pressure,
                                                                const double       delta_forward_gibbs_energy,
                                                                const double       cumulative_forward_reaction_progress,
                                                                const unsigned int reaction_index) const
      {
        AssertIndexRange(reaction_index, kinetic_prefactors.size());

        AssertThrow(temperature > 0.0,
                    ExcMessage("Temperature must be strictly positive in net_forward_reaction_rate(), "
                               "but received " +
                               std::to_string(temperature) + " K."));

        const double activation_enthalpy = activation_enthalpies[reaction_index];
        const double activation_volume   = activation_volumes[reaction_index];
        AssertThrow(activation_enthalpy >= 0.0,
                    ExcMessage("Activation enthalpy must be non-negative, "
                               "but received " +
                               std::to_string(activation_enthalpy) + " J/mol."));

        const double arrhenius_factor =
          std::exp(-(activation_enthalpy + (pressure * activation_volume)) / (constants::gas_constant * temperature));
        const double magnitude            = 1.0 - std::exp(-std::abs(delta_forward_gibbs_energy) / (constants::gas_constant * temperature));
        const double thermodynamic_factor = -std::copysign(1.0, delta_forward_gibbs_energy) * magnitude;
        const double reaction_progress =
          (delta_forward_gibbs_energy < 0.0) ? (1.0 - cumulative_forward_reaction_progress) : cumulative_forward_reaction_progress;

        return kinetic_prefactors[reaction_index] * temperature * arrhenius_factor * thermodynamic_factor * reaction_progress;
      }

      template <int dim>
      void
      InterfaceControlledGrowth<dim>::declare_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Interface controlled growth");
        {
          prm.declare_entry("Kinetic prefactors",
                            "7.0e7",
                            Patterns::List(Patterns::Double(0.0), 1, Patterns::List::max_int_value, "|"),
                            "Kinetic prefactors Z, one '|'-separated entry per reaction assigned to this model. Units: 1/s/K");
          prm.declare_entry("Activation enthalpies",
                            "274e3",
                            Patterns::List(Patterns::Double(0.0), 1, Patterns::List::max_int_value, "|"),
                            "Activation enthalpies Ha, one '|'-separated entry per reaction assigned to this model. Units: J/mol");
          prm.declare_entry("Activation volumes",
                            "3.3e-6",
                            Patterns::List(Patterns::Double(0.0), 1, Patterns::List::max_int_value, "|"),
                            "Activation volumes Va, one '|'-separated entry per reaction assigned to this model. Units: m^3/mol");
        }
        prm.leave_subsection();
      }

      template <int dim>
      void
      InterfaceControlledGrowth<dim>::parse_parameters(ParameterHandler &prm, const unsigned int n_reactions)
      {
        prm.enter_subsection("Interface controlled growth");
        {
          kinetic_prefactors    = Utilities::string_to_double(Utilities::split_string_list(prm.get("Kinetic prefactors"), '|'));
          activation_enthalpies = Utilities::string_to_double(Utilities::split_string_list(prm.get("Activation enthalpies"), '|'));
          activation_volumes    = Utilities::string_to_double(Utilities::split_string_list(prm.get("Activation volumes"), '|'));

          AssertThrow(kinetic_prefactors.size() == n_reactions && activation_enthalpies.size() == n_reactions &&
                        activation_volumes.size() == n_reactions,
                      ExcMessage("'Interface controlled growth' parameter lists must each have exactly " + std::to_string(n_reactions) +
                                 " '|'-separated entries, matching the number of "
                                 "reactions assigned to this kinetics model."));
        }
        prm.leave_subsection();
      }
    }
  }
}

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      ASPECT_REGISTER_REACTION_MODEL(
        InterfaceControlledGrowth,
        "Interface controlled growth",
        "Interface-controlled growth kinetics following Cahn (1956; https://doi.org/10.1016/0001-6160(56)90041-4) and "
        "Hosoya et al. (2005; https://doi.org/10.1029/2005GL023398)")
    }
  }
}
