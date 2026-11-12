/*
  Copyright (C) 2022 - 2026 by the authors of the ASPECT code.

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

#include <aspect/simulator_signals.h>
#include <aspect/simulator_access.h>

using namespace aspect;


template <int dim>
void modify_reactions (const SimulatorAccess<dim> &simulator_access,
                       const MaterialModel::MaterialModelInputs<dim> &in,
                       MaterialModel::MaterialModelOutputs<dim> &out)
{
  // If the current cell being processed in compute_reactions is not valid or locally owned, skip.
  if (in.current_cell.state() != IteratorState::valid || !in.current_cell->is_locally_owned())
    return;

  // Get the reaction rate outputs so that we can modify it's entries
  const std::shared_ptr<MaterialModel::ReactionRateOutputs<dim>> reaction_rate_outputs
    = out.template get_additional_output_object<MaterialModel::ReactionRateOutputs<dim>>();
  const unsigned int porosity_idx = simulator_access.introspection().compositional_index_for_name("porosity");
  const unsigned int bound_fluid_idx = simulator_access.introspection().compositional_index_for_name("bound_fluid");

  for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
    {
      reaction_rate_outputs->reaction_rates[q][porosity_idx] = 0;
      reaction_rate_outputs->reaction_rates[q][bound_fluid_idx] = 0;
    }
}


template <int dim>
void signal_connector (SimulatorSignals<dim> &signals)
{
  signals.modify_reactions.connect (&modify_reactions<dim>);
}


ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                  signal_connector<3>)
