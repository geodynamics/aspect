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

#include <aspect/material_model/reaction_model/reaction_chain.h>
#include <aspect/global.h>
#include <algorithm>
#include <map>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      template <int dim>
      void
      ReactionChain<dim>::initialize()
      {
        for (auto &kinetics : kinetics_models)
          {
            kinetics->initialize_simulator(this->get_simulator());
          }
      }

      template <int dim>
      unsigned int
      ReactionChain<dim>::n_reactions() const
      {
        return reactions.size();
      }

      template <int dim>
      std::vector<double>
      ReactionChain<dim>::clamp_cumulative_progress(std::vector<double> reaction_progress_values) const
      {
        AssertThrow(reaction_progress_values.size() == reactions.size(),
                    ExcMessage("Size of 'reaction_progress_values' (" + std::to_string(reaction_progress_values.size()) +
                               ") does not match number of reactions (" + std::to_string(reactions.size()) + ")."));

        // Sanity check: reaction progress values should be close to [0, 1]
        for (unsigned int i = 0; i < reaction_progress_values.size(); ++i)
          {
            const double val = reaction_progress_values[i];

            AssertThrow(std::isfinite(val),
                        ExcMessage("'reaction_progress_values[" + std::to_string(i) + "]' = " + std::to_string(val) +
                                   " is not finite (NaN or Inf)."));

            AssertThrow(val > -tolerance_in_reaction_progress && val < 1.0 + tolerance_in_reaction_progress,
                        ExcMessage("'reaction_progress_values[" + std::to_string(i) + "]' = " + std::to_string(val) +
                                   " is too far outside the physical range [0, 1] (allowed excursion is " +
                                   std::to_string(tolerance_in_reaction_progress) + ")."));
          }

        // Sanity check: reaction progress values should be close to monotonic
        for (unsigned int i = 1; i < reaction_progress_values.size(); ++i)
          AssertThrow(reaction_progress_values[i] < reaction_progress_values[i-1] + tolerance_in_reaction_progress,
                      ExcMessage("'reaction_progress_values[" + std::to_string(i) + "]' = " +
                                 std::to_string(reaction_progress_values[i]) + " exceeds "
                                 "'reaction_progress_values[" + std::to_string(i-1) + "]' = " +
                                 std::to_string(reaction_progress_values[i-1]) + " by more than the "
                                 "allowed excursion (" + std::to_string(tolerance_in_reaction_progress) + "). Reaction " +
                                 std::to_string(i) + " should not have progressed further than reaction " +
                                 std::to_string(i-1) + " in a sequential chain."));

        if (!reaction_progress_values.empty())
          reaction_progress_values[0] = std::clamp(reaction_progress_values[0], 0.0, 1.0);

        for (unsigned int i = 1; i < reaction_progress_values.size(); ++i)
          {
            reaction_progress_values[i] = std::clamp(reaction_progress_values[i], 0.0, 1.0);  // Clamp to valid physical range [0, 1]
            reaction_progress_values[i] = std::min(reaction_progress_values[i], reaction_progress_values[i-1]);  // Enforce monotonicity xi[i] <= xi[i-1]
          }

        return reaction_progress_values;
      }

      template <int dim>
      std::vector<double>
      ReactionChain<dim>::compute_phase_mass_fractions(const std::vector<double> &reaction_progress_values, const std::vector<double> &phase_densities) const
      {
        AssertThrow(!reaction_progress_values.empty(),
                    ExcMessage("'reaction_progress_values' cannot be empty."));
        AssertThrow(reaction_progress_values.size() == reactions.size(),
                    ExcMessage("Size of 'reaction_progress_values' (" + std::to_string(reaction_progress_values.size()) +
                               ") does not match number of reactions (" + std::to_string(reactions.size()) + ")."));
        AssertThrow(phase_densities.size() == reactions.size() + 1,
                    ExcMessage("Size of 'phase_densities' (" + std::to_string(phase_densities.size()) +
                               ") must be equal to reactions + 1 (" + std::to_string(reactions.size() + 1) + ")."));

        for (unsigned int p = 0; p < phase_densities.size(); ++p)
          AssertThrow(phase_densities[p] > 0.0 && std::isfinite(phase_densities[p]),
                      ExcMessage("Phase density at index " + std::to_string(p) + " must be strictly positive and finite."));

        std::vector<double> volume_fraction_values(phase_densities.size());
        std::vector<double> mass_fraction_values(phase_densities.size(), 0.0);

        // Guarantee bounds and monotonicity of reaction sequence
        const std::vector<double> reaction_progress_values_clamped = clamp_cumulative_progress(reaction_progress_values);

        // Derive (N+1) volume fractions from (N) cumulative reaction progress fields
        volume_fraction_values[0] = 1.0 - reaction_progress_values_clamped[0];
        for (unsigned int i = 1; i < reaction_progress_values_clamped.size(); ++i)
          volume_fraction_values[i] = std::max(0.0, reaction_progress_values_clamped[i-1] - reaction_progress_values_clamped[i]);
        volume_fraction_values.back() = std::max(0.0, reaction_progress_values_clamped.back());

        // Sanity check: by construction these should sum to 1 and be non-negative
        const double volume_fraction_sum = std::accumulate(volume_fraction_values.begin(), volume_fraction_values.end(), 0.0);
        AssertThrow(std::abs(volume_fraction_sum - 1.0) < 1e-10,
                    ExcMessage("Computed phase volume fractions sum to " + std::to_string(volume_fraction_sum) + " instead of 1."));

        // Compute bulk density of phase mixture
        double bulk_density = 0.0;
        for (unsigned int p = 0; p < phase_densities.size(); ++p)
          bulk_density += phase_densities[p] * volume_fraction_values[p];

        // Convert (N+1) phase volume fractions to (N+1) phase mass fractions
        if (bulk_density > 0.0)
          for (unsigned int p = 0; p < phase_densities.size(); ++p)
            mass_fraction_values[p] = (phase_densities[p] * volume_fraction_values[p]) / bulk_density;

        return mass_fraction_values;
      }

      template <int dim>
      double
      ReactionChain<dim>::net_forward_reaction_rate(const double temperature,
                                                    const double pressure,
                                                    const double delta_forward_gibbs_energy,
                                                    const double cumulative_forward_reaction_progress,
                                                    const unsigned int reaction_index) const
      {
        AssertIndexRange(reaction_index, reactions.size());
        const ReactionStep<dim> &reaction = reactions[reaction_index];
        return reaction.kinetics->net_forward_reaction_rate(temperature, pressure,
                                                            delta_forward_gibbs_energy,
                                                            cumulative_forward_reaction_progress,
                                                            reaction.local_reaction_index);
      }

      template <int dim>
      void
      ReactionChain<dim>::declare_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Reaction chain");
        {
          prm.declare_entry("Kinetic models",
                            "Interface controlled growth",
                            Patterns::List(Patterns::Selection(ReactionModelPluginList<dim>::get_pattern_of_names()),
                                           1, Patterns::List::max_int_value, "|"),
                            "'|'-separated list of registered reaction kinetics model names, one per reaction in the chain. The number of entries "
                            "sets the number of reactions N, which must equal the number of compositional fields tracking reaction progress.");

          prm.declare_entry("Tolerance in reaction progress",
                            "0.1",
                            Patterns::Double(0.0),
                            "Allowed excursion before a reaction progress value outside [0, 1] is considered unphysical.");

          // One subsection per registered kinetics model (not per reaction).
          ReactionModelPluginList<dim>::declare_parameters(prm);
        }
        prm.leave_subsection();
      }

      template <int dim>
      void
      ReactionChain<dim>::parse_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Reaction chain");
        {
          tolerance_in_reaction_progress = prm.get_double("Tolerance in reaction progress");

          AssertThrow(tolerance_in_reaction_progress >= 0.0 && tolerance_in_reaction_progress < 1.0,
                      ExcMessage("The 'Tolerance in reaction progress' must be between [0, 1), "
                                 "but was set to " + std::to_string(tolerance_in_reaction_progress) + "."));

          const std::vector<std::string> kinetic_model_names = Utilities::split_string_list(prm.get("Kinetic models"), '|');
          const unsigned int n_reactions = kinetic_model_names.size();

          AssertThrow(n_reactions > 0, ExcMessage("'Reaction chain/Kinetic models' must contain at least one entry."));

          reactions.resize(n_reactions);
          kinetics_models.clear();

          // Group global reaction indices by model name, preserving first-appearance order, so reactions sharing a model name get consecutive
          // local indices matching the order that name appears in Kinetic models.
          std::vector<std::string> unique_model_names;
          std::map<std::string, std::vector<unsigned int>> global_indices_by_model;
          for (unsigned int global_index = 0; global_index < n_reactions; ++global_index)
            {
              const std::string &model_name = kinetic_model_names[global_index];
              if (global_indices_by_model.find(model_name) == global_indices_by_model.end())
                unique_model_names.push_back(model_name);
              global_indices_by_model[model_name].push_back(global_index);
            }

          for (const std::string &model_name : unique_model_names)
            {
              const std::vector<unsigned int> &global_indices = global_indices_by_model[model_name];

              std::shared_ptr<Cahn1956Interface<dim>> kinetics(create_reaction_model<dim>(model_name).release());
              kinetics->parse_parameters(prm, global_indices.size());
              kinetics_models.push_back(kinetics);

              for (unsigned int local_index = 0; local_index < global_indices.size(); ++local_index)
                {
                  reactions[global_indices[local_index]].kinetics             = kinetics;
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
#define INSTANTIATE(dim) template struct ReactionStep<dim>; template class ReactionChain<dim>;
      ASPECT_INSTANTIATE(INSTANTIATE)
#undef INSTANTIATE
    }
  }
}
