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

#include <aspect/material_model/reaction_model/gibbs_free_energy/linearized.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      namespace GibbsFreeEnergyModel
      {
        template <int dim>
        double Linearized<dim>::compute_delta_gibbs_free_energy(const double temperature,
                                                                const double pressure,
                                                                const GibbsFreeEnergyInputs &reference_state,
                                                                const unsigned int) const
        {
          return reference_state.reference_delta_gibbs_free_energy
                 + (pressure - reference_state.reference_pressure) * reference_state.delta_volume
                 - (temperature - reference_state.reference_temperature) * reference_state.delta_entropy;
        }



        template <int dim>
        void Linearized<dim>::declare_parameters(ParameterHandler &)
        {}



        template <int dim>
        void Linearized<dim>::parse_parameters(ParameterHandler &, const unsigned int)
        {}
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
      namespace GibbsFreeEnergyModel
      {
        ASPECT_REGISTER_GIBBS_FREE_ENERGY_MODEL(Linearized,
                                                "Linearized Gibbs free energy",
                                                "Computes G_B-G_A from a reference value using linear pressure-volume and temperature-entropy corrections.")
      }
    }
  }
}
