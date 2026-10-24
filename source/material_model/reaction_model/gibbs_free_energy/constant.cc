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

#include <aspect/material_model/reaction_model/gibbs_free_energy/constant.h>

#include <deal.II/base/patterns.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      namespace GibbsFreeEnergyModel
      {
        template <int dim>
        double Constant<dim>::compute_delta_gibbs_free_energy(const double,
                                                              const double,
                                                              const GibbsFreeEnergyInputs &,
                                                              const unsigned int reaction_index) const
        {
          AssertIndexRange(reaction_index, delta_gibbs_free_energies.size());
          return delta_gibbs_free_energies[reaction_index];
        }



        template <int dim>
        void Constant<dim>::declare_parameters(ParameterHandler &prm)
        {
          prm.enter_subsection("Constant Gibbs free energy");
          {
            prm.declare_entry("Gibbs free energy differences",
                              "0.0",
                              Patterns::List(Patterns::Double(), 1, Patterns::List::max_int_value, "|"),
                              "Constant values of G_B-G_A, one '|'-separated entry per reaction assigned to this model. "
                              "Negative values favor A -> B and positive values favor B -> A. Units: J/mol.");
          }
          prm.leave_subsection();
        }



        template <int dim>
        void Constant<dim>::parse_parameters(ParameterHandler &prm, const unsigned int n_reactions)
        {
          prm.enter_subsection("Constant Gibbs free energy");
          {
            delta_gibbs_free_energies =
              Utilities::string_to_double(Utilities::split_string_list(prm.get("Gibbs free energy differences"), '|'));

            AssertThrow(delta_gibbs_free_energies.size() == n_reactions,
                        ExcMessage("The 'Constant Gibbs free energy/Gibbs free energy differences' parameter must have exactly " +
                                   std::to_string(n_reactions) + " '|'-separated entries, matching the number of reactions assigned to this model."));
          }
          prm.leave_subsection();
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
      namespace GibbsFreeEnergyModel
      {
        ASPECT_REGISTER_GIBBS_FREE_ENERGY_MODEL(Constant,
                                                "Constant Gibbs free energy",
                                                "Returns a constant Gibbs free-energy difference G_B-G_A for each reaction.")
      }
    }
  }
}
