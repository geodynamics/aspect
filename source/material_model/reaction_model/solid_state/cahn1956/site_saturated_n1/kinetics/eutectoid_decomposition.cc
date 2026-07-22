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

#include <aspect/material_model/reaction_model/solid_state/cahn1956/site_saturated_n1/kinetics/eutectoid_decomposition.h>
#include <aspect/global.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/patterns.h>

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
            double EutectoidDecomposition<dim>::clamped_exp(const double x)
            {
              return std::exp(std::clamp(x, -700.0, 700.0));
            }



            template <int dim>
            double EutectoidDecomposition<dim>::arrhenius_factor(const double temperature, const double /*pressure*/) const
            {
              const double safe_temperature = std::max(temperature, 1.0e-10);
              return clamped_exp(-activation_energy / (gas_constant * safe_temperature));
            }



            template <int dim>
            double EutectoidDecomposition<dim>::thermodynamic_factor(const double /*temperature*/, const double delta_gibbs_energy) const
            {
              return -delta_gibbs_energy * std::abs(delta_gibbs_energy);
            }



            template <int dim>
            double EutectoidDecomposition<dim>::
            reaction_rate(const double temperature, const double pressure, const double delta_gibbs_energy, const double reactant_phase_fraction) const
            {
              return kinetic_factor * thermodynamic_factor(temperature, delta_gibbs_energy) * arrhenius_factor(temperature, pressure) * reactant_phase_fraction;
            }



            template <int dim>
            void EutectoidDecomposition<dim>::declare_parameters(ParameterHandler &prm)
            {
              prm.enter_subsection("Eutectoid decomposition");
              {
                prm.declare_entry("Kinetic factor",
                                  "2.7e-16",
                                  Patterns::Double(0.0),
                                  "Kinetic prefactor Z in the eutectoid decomposition lamellae growth law "
                                  "dX/dt = Z * (-dG) * |dG| * exp(-Ea/(R*T)) * X. "
                                  "Units: mol^2/J^2/s");
                prm.declare_entry("Activation energy",
                                  "355e3",
                                  Patterns::Double(0.0),
                                  "Activation energy Ea of diffusion. Units: J/mol");
              }
              prm.leave_subsection();
            }



            template <int dim>
            void EutectoidDecomposition<dim>::parse_parameters(ParameterHandler &prm)
            {
              prm.enter_subsection("Eutectoid decomposition");
              {
                kinetic_factor    = prm.get_double("Kinetic factor");
                activation_energy = prm.get_double("Activation energy");
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
              EutectoidDecomposition,
              "Eutectoid decomposition",
              "Eutectoid decomposition kinetics following Kubo et al. (2000, 2002): "
              "dX/dt = Z * (-dG) * |dG| * exp(-Ea/(R*T)) * X.")
          }
        }
      }
    }
  }
}
