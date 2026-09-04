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

#include <aspect/material_model/reaction_model/kinetics/constant.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/patterns.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      template <int dim>
      double ConstantReactionRate<dim>::
      net_forward_reaction_rate(const double,
                                const double,
                                const double,
                                const double,
                                const unsigned int reaction_index) const
      {
        AssertIndexRange(reaction_index, reaction_rates.size());
        return reaction_rates[reaction_index];
      }

      template <int dim>
      void ConstantReactionRate<dim>::declare_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Constant reaction rate");
        {
          prm.declare_entry("Reaction rates",
                            "0.0",
                            Patterns::List(Patterns::Double(), 1, Patterns::List::max_int_value, ","),
                            "Constant net forward reaction rates, one comma-separated entry per reaction assigned to this model. "
                            "Units: 1/s, or 1/yr if the global 'Use years instead of seconds' parameter is set.");
        }
        prm.leave_subsection();
      }

      template <int dim>
      void ConstantReactionRate<dim>::parse_parameters(ParameterHandler &prm, const unsigned int n_reactions)
      {
        prm.enter_subsection("Constant reaction rate");
        {
          reaction_rates = Utilities::string_to_double(Utilities::split_string_list(prm.get("Reaction rates"), ','));

          AssertThrow(reaction_rates.size() == n_reactions,
                      ExcMessage("The 'Constant reaction rate/Reaction rates' parameter must have exactly " +
                                 std::to_string(n_reactions) + " comma-separated entries, matching the number of "
                                 "reactions assigned to this kinetics model."));

          if (this->convert_output_to_years())
            for (double &reaction_rate : reaction_rates)
              reaction_rate /= year_in_seconds;
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
        ConstantReactionRate,
        "Constant reaction rate",
        "A trivial reaction kinetics model that returns a user-specified constant net forward reaction rate for testing purposes.")
    }
  }
}
