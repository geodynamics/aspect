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

#include <aspect/material_model/reaction_model/solid_state/reaction_chain.h>
#include <aspect/global.h>
#include <algorithm>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      namespace SolidState
      {
        template <int dim>
        bool ReactionChain<dim>::any_enabled() const
        {
          return std::any_of(reactions.begin(), reactions.end(), [](const ReactionStep<dim> &r)
          {
            return r.enabled;
          });
        }

        template <int dim>
        std::vector<double>
        ReactionChain<dim>::clamp_cumulative_progress(std::vector<double> xi) const
        {
          Assert(xi.size() == reactions.size(), ExcInternalError());

          if (!xi.empty())
            xi[0] = std::clamp(xi[0], 0.0, 1.0);

          for (unsigned int i = 1; i < xi.size(); ++i)
            xi[i] = std::clamp(xi[i], 0.0, xi[i-1]);

          return xi;
        }

        template <int dim>
        std::vector<double>
        ReactionChain<dim>::compute_phase_mass_fractions(const std::vector<double> &xi_raw) const
        {
          Assert(xi_raw.size() + 1 == phase_names.size(), ExcInternalError());

          if (xi_raw.empty())
            return {1.0};

          // Guarantee bounds and monotonicity
          const std::vector<double> xi = clamp_cumulative_progress(xi_raw);

          std::vector<double> X(phase_names.size());
          X.front() = std::max(0.0, 1.0 - xi.front());
          for (unsigned int i = 1; i < xi.size(); ++i)
            X[i] = std::max(0.0, xi[i-1] - xi[i]);
          X.back() = std::max(0.0, xi.back());
          return X;
        }

        template <int dim>
        void
        ReactionChain<dim>::declare_parameters(ParameterHandler &prm)
        {
          prm.enter_subsection("Reaction chain");
          {
            const std::string default_phases = "olivine, wadsleyite, ringwoodite, postspinel";
            prm.declare_entry("Phase names",
                              default_phases,
                              Patterns::List(Patterns::Anything()),
                              "Comma-separated list of phase names defining the sequential chain.");
            const std::vector<std::string> default_phase_names = Utilities::split_string_list(default_phases);

            for (unsigned int i = 0; i + 1 < default_phase_names.size(); ++i)
              {
                prm.enter_subsection(default_phase_names[i] + " " + default_phase_names[i+1]);
                {
                  prm.declare_entry("Enable transformation", "true", Patterns::Bool(),
                                    "Whether to enable the " + default_phase_names[i] + " <-> " + default_phase_names[i+1] + " transformation.");
                  Utilities::AsciiDataProfile<dim>::declare_parameters(prm, "$ASPECT_SOURCE_DIR/data/material-model/reaction-chain/", "profile.tsv");
                  prm.declare_entry("Data directory", "$ASPECT_SOURCE_DIR/data/material-model/reaction-chain/",
                                    Patterns::DirectoryName(), "Directory containing this transformation's thermodynamic data file.");
                  prm.declare_entry("Data file name", "profile.tsv", Patterns::Anything(),
                                    "Name of this transformation's thermodynamic data file.");
                  prm.declare_entry("Kinetics model", "Interface controlled growth", Patterns::Anything(),
                                    "Registered Cahn1956::SiteSaturatedN1 reaction kinetics plugin governing this transformation.");
                  Cahn1956::SiteSaturatedN1::declare_parameters<dim>(prm);
                }
                prm.leave_subsection();
              }
          }
          prm.leave_subsection();
        }

        template <int dim>
        void
        ReactionChain<dim>::parse_parameters(ParameterHandler &prm)
        {
          phase_names = Utilities::split_string_list(prm.get("Phase names"));
          reactions.resize(phase_names.size() - 1);

          for (unsigned int i = 0; i < reactions.size(); ++i)
            {
              prm.enter_subsection(phase_names[i] + " " + phase_names[i+1]);
              {
                auto &step = reactions[i];
                step.enabled = prm.get_bool("Enable transformation");
                step.cumulative_field_name = "xi_" + phase_names[i] + "_" + phase_names[i+1];

                if (step.enabled)
                  {
                    step.data_directory = Utilities::expand_ASPECT_SOURCE_DIR(prm.get("Data directory"));
                    step.data_filename  = prm.get("Data file name");
                    prm.enter_subsection("Ascii data model");
                    {
                      prm.set("Data directory", step.data_directory);
                      prm.set("Data file name", step.data_filename);
                    }
                    prm.leave_subsection();
                    step.profile.parse_parameters(prm);

                    step.kinetics = Cahn1956::SiteSaturatedN1::create_reaction_model<dim>(prm.get("Kinetics model"));
                    step.kinetics->parse_parameters(prm);
                  }
              }
              prm.leave_subsection();
            }
        }

        template <int dim>
        void
        ReactionChain<dim>::initialize()
        {
          for (auto &step : reactions)
            if (step.enabled)
              {
                step.profile.initialize(this->get_mpi_communicator());
                step.kinetics->initialize_simulator(this->get_simulator());
                step.dG_idx = step.profile.get_column_index_from_name("delta_molar_gibbs");
                step.dS_idx = step.profile.get_column_index_from_name("delta_molar_entropy");
                step.dV_idx = step.profile.get_column_index_from_name("delta_molar_volume");
              }
        }

        template <int dim>
        const std::vector<std::string> &ReactionChain<dim>::get_phase_names() const
        {
          return phase_names;
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
#define INSTANTIATE(dim) template struct ReactionStep<dim>; template class ReactionChain<dim>;
        ASPECT_INSTANTIATE(INSTANTIATE)
#undef INSTANTIATE
      }
    }
  }
}
