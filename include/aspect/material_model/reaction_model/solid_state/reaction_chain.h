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

#ifndef _aspect_material_model_reaction_model_solid_state_reaction_chain_h
#define _aspect_material_model_reaction_model_solid_state_reaction_chain_h

#include <aspect/material_model/reaction_model/solid_state/cahn1956/site_saturated_n1/interface.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      namespace SolidState
      {
        template <int dim>
        struct ReactionStep
        {
          // Thermodynamic data (adiabatic profile)
          Utilities::AsciiDataProfile<dim> profile;

          // Kinetic model.
          //
          // NOTE: this currently hardcodes the Cahn1956::SiteSaturatedN1 kinetics family. ReactionChain does not yet support mixing kinetics
          // families with different Interface signatures (e.g. a future JMAK-type or general Cahn1956 model) within a single chain. When
          // such a family is added, this will need to become either a type-erased/variant kinetics holder, or ReactionChain (and ReactionStep)
          // will need to be templated on the kinetics family's Interface type.
          std::unique_ptr<Cahn1956::SiteSaturatedN1::Interface<dim>> kinetics;

          // Thermodynamic data file
          std::string data_directory;
          std::string data_filename;

          // Cumulative reaction progress composition field (e.g., "xi_olivine_wadsleyite")
          std::string cumulative_field_name;

          // Column indices of thermodynamic data
          unsigned int dG_idx = numbers::invalid_unsigned_int;
          unsigned int dS_idx = numbers::invalid_unsigned_int;
          unsigned int dV_idx = numbers::invalid_unsigned_int;

          // Enable reaction
          bool enabled = true;
        };

        /**
         * A linear chain of N polymorphic transformations among N+1 phases, driven by N cumulative reaction-progress compositional fields with
         * enforced ordering 1 >= xi_0 >= ... >= xi_{N-1} >= 0.
         *
         * Scope: independent linear chains only. Branching networks (a phase produced/consumed by more than one reaction) or competing back-
         * reaction pathways need a stoichiometry-matrix mass balance and are out of scope here.
         */
        template <int dim>
        class ReactionChain : public ::aspect::SimulatorAccess<dim>
        {
          public:
            std::vector<std::string> phase_names;
            std::vector<ReactionStep<dim>> reactions;

            bool any_enabled() const;

            std::vector<double> clamp_cumulative_progress(std::vector<double> xi_cumulative) const;

            std::vector<double> compute_phase_mass_fractions(const std::vector<double> &xi_clamped) const;

            static void declare_parameters(ParameterHandler &prm);

            void parse_parameters(ParameterHandler &prm);

            void initialize();

            const std::vector<std::string> &get_phase_names() const;
        };
      }
    }
  }
}

#endif
