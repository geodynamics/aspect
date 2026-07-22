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

#include <aspect/material_model/reaction_model/solid_state/cahn1956/site_saturated_n1/kinetics/interface_controlled_growth.h>
#include <aspect/global.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/patterns.h>

#include <algorithm>
#include <cmath>

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
            template <int dim>
            double InterfaceControlledGrowth<dim>::clamped_exp(const double x)
            {
              return std::exp(std::clamp(x, -700.0, 700.0));
            }



            template <int dim>
            double InterfaceControlledGrowth<dim>::arrhenius_factor(const double temperature, const double pressure) const
            {
              const double safe_temperature = std::max(temperature, 1.0e-10);
              return clamped_exp(-(activation_enthalpy + (pressure * activation_volume)) / (gas_constant * safe_temperature));
            }



            template <int dim>
            double InterfaceControlledGrowth<dim>::thermodynamic_factor(const double temperature, const double delta_gibbs_energy) const
            {
              const double safe_temperature = std::max(temperature, 1.0e-10);
              const double magnitude = 1.0 - clamped_exp(-std::abs(delta_gibbs_energy) / (gas_constant * safe_temperature));
              return -std::copysign(1.0, delta_gibbs_energy) * magnitude;
            }



            template <int dim>
            double InterfaceControlledGrowth<dim>::
            reaction_rate(const double temperature, const double pressure, const double delta_gibbs_energy, const double reactant_phase_fraction) const
            {
              return kinetic_factor * temperature * arrhenius_factor(temperature, pressure) * thermodynamic_factor(temperature, delta_gibbs_energy) * reactant_phase_fraction;
            }



            template <int dim>
            void InterfaceControlledGrowth<dim>::declare_parameters(ParameterHandler &prm)
            {
              prm.enter_subsection("Interface controlled growth");
              {
                prm.declare_entry("Kinetic factor",
                                  "7.0e7",
                                  Patterns::Double(0.0),
                                  "Kinetic prefactor Z in the interface-controlled growth law "
                                  "dX/dt = Z * T * exp(-(Ha + P*Va)/(R*T)) * (1 - exp(-|dG|/(R*T))) * X. "
                                  "Units: 1/s/K");
                prm.declare_entry("Activation enthalpy",
                                  "274e3",
                                  Patterns::Double(0.0),
                                  "Activation enthalpy Ha of the interface reaction. Units: J/mol");
                prm.declare_entry("Activation volume",
                                  "3.3e-6",
                                  Patterns::Double(0.0),
                                  "Activation volume Va of the interface reaction. Units: m^3/mol");
              }
              prm.leave_subsection();
            }



            template <int dim>
            void InterfaceControlledGrowth<dim>::parse_parameters(ParameterHandler &prm)
            {
              prm.enter_subsection("Interface controlled growth");
              {
                kinetic_factor       = prm.get_double("Kinetic factor");
                activation_enthalpy  = prm.get_double("Activation enthalpy");
                activation_volume    = prm.get_double("Activation volume");
              }
              prm.leave_subsection();
            }
          }
        }
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
      namespace SolidState
      {
        namespace Cahn1956
        {
          namespace SiteSaturatedN1
          {
            ASPECT_REGISTER_CAHN1956_SITE_SATURATED_N1_REACTION_MODEL(
              InterfaceControlledGrowth,
              "Interface controlled growth",
              "Interface-controlled polymorphic transformation kinetics following Hosoya et al. (2005): "
              "dX/dt = Z * T * exp(-(Ha + P*Va)/(R*T)) * (1 - exp(-|dG|/(R*T))) * X.")
          }
        }
      }
    }
  }
}
