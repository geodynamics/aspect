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

#ifndef _aspect_material_model_reaction_model_reaction_chain_h
#define _aspect_material_model_reaction_model_reaction_chain_h

#include <aspect/material_model/reaction_model/kinetics/cahn1956_interface.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      /**
       * Stores the kinetics model for each reaction in a ReactionChain object. Each distinct reaction is referenced by a local index.
       */
      template <int dim>
      struct ReactionStep
      {
        std::shared_ptr<Cahn1956Interface<dim>> kinetics;
        unsigned int local_reaction_index = numbers::invalid_unsigned_int;
      };

      /**
       * A linear chain of N phase transformations among N+1 phases, driven by N cumulative reaction progress fields with enforced ordering
       *
       *   1 >= xi_0 >= ... >= xi_{N-1} >= 0.
       *
       * N is set by the number of entries in "Reaction chain/Kinetic models" (see Cahn1956Interface).
       *
       * Independent linear chains only. Branching networks (a phase produced/consumed by more than one reaction) or competing back-reaction
       * pathways need a stoichiometry-matrix mass balance and are out of scope here.
       */
      template <int dim>
      class ReactionChain : public ::aspect::SimulatorAccess<dim>
      {
        public:
          std::vector<ReactionStep<dim>> reactions;

          void initialize();
          unsigned int n_reactions() const;

          std::vector<double> clamp_cumulative_progress(std::vector<double> reaction_progress_values) const;
          std::vector<double> compute_phase_mass_fractions(const std::vector<double> &reaction_progress_values,
                                                           const std::vector<double> &phase_densities) const;

          double net_forward_reaction_rate(const double temperature,
                                           const double pressure,
                                           const double delta_forward_gibbs_energy,
                                           const double cumulative_forward_reaction_progress,
                                           const unsigned int reaction_index) const;

          static void declare_parameters(ParameterHandler &prm);
          void parse_parameters(ParameterHandler &prm);

        private:
          /**
           * One entry per unique kinetics model name used in the chain.
           */
          std::vector<std::shared_ptr<Cahn1956Interface<dim>>> kinetics_models;

          /**
           * Allowed excursion/slack before a reaction progress value is considered too extreme for a given timestep.
           */
          double tolerance_in_reaction_progress;
      };
    }
  }
}
#endif
