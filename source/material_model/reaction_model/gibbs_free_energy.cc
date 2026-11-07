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

#include <aspect/material_model/reaction_model/gibbs_free_energy.h>
#include <aspect/global.h>

#include <map>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      template <int dim>
      double GibbsFreeEnergy<dim>::delta_gibbs_free_energy(const double temperature,
                                                           const double pressure,
                                                           const GibbsFreeEnergyInputs &reference_state,
                                                           const unsigned int reaction_index) const
      {
        AssertIndexRange(reaction_index, reactions.size());

        const GibbsFreeEnergyStep<dim> &reaction = reactions[reaction_index];
        const double delta_gibbs = reaction.model->compute_delta_gibbs_free_energy(temperature,
                                                                                   pressure,
                                                                                   reference_state,
                                                                                   reaction.local_reaction_index);

        AssertThrow(std::isfinite(delta_gibbs), ExcMessage("The computed Gibbs free-energy difference must be finite."));

        return delta_gibbs;
      }



      template <int dim>
      unsigned int GibbsFreeEnergy<dim>::n_reactions() const
      {
        return reactions.size();
      }



      template <int dim>
      void GibbsFreeEnergy<dim>::declare_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Gibbs free energy");
        {
          prm.declare_entry("Models",
                            "Constant Gibbs free energy",
                            Patterns::List(Patterns::Anything()),
                            "A '|'-separated Gibbs free-energy model for each reaction in the sequential chain. "
                            "Each model computes G_B-G_A for adjacent phases A and B.");
          GibbsFreeEnergyModel::PluginList<dim>::declare_parameters(prm);
        }
        prm.leave_subsection();
      }



      template <int dim>
      void GibbsFreeEnergy<dim>::parse_parameters(ParameterHandler &prm, const unsigned int n_reactions)
      {
        prm.enter_subsection("Gibbs free energy");
        {
          const std::vector<std::string> model_names = Utilities::split_string_list(prm.get("Models"), '|');
          AssertThrow(model_names.size() == n_reactions,
                      ExcMessage("The number of entries in 'Gibbs free energy/Models' must match the number of reactions."));

          reactions.resize(n_reactions);

          std::vector<std::string> unique_model_names;
          std::map<std::string, std::vector<unsigned int>> global_indices_by_model;
          for (unsigned int global_index = 0; global_index < n_reactions; ++global_index)
            {
              const std::string &model_name = model_names[global_index];
              if (global_indices_by_model.find(model_name) == global_indices_by_model.end())
                unique_model_names.push_back(model_name);
              global_indices_by_model[model_name].push_back(global_index);
            }

          for (const std::string &model_name : unique_model_names)
            {
              const std::vector<unsigned int> &global_indices = global_indices_by_model[model_name];
              std::shared_ptr<GibbsFreeEnergyModel::Interface<dim>> model(
                GibbsFreeEnergyModel::create_gibbs_free_energy_model<dim>(model_name).release());
              model->parse_parameters(prm, global_indices.size());

              for (unsigned int local_index = 0; local_index < global_indices.size(); ++local_index)
                {
                  reactions[global_indices[local_index]].model = model;
                  reactions[global_indices[local_index]].local_reaction_index = local_index;
                }
            }
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
#define INSTANTIATE(dim) template class GibbsFreeEnergy<dim>;
      ASPECT_INSTANTIATE(INSTANTIATE)
#undef INSTANTIATE
    }
  }
}
