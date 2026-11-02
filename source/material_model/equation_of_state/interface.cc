// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/material_model/equation_of_state/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    EquationOfStateOutputs<dim>::EquationOfStateOutputs(const unsigned int n_individual_compositions_and_phases)
      :
      densities(n_individual_compositions_and_phases, numbers::signaling_nan<double>()),
      thermal_expansion_coefficients(n_individual_compositions_and_phases, numbers::signaling_nan<double>()),
      specific_heat_capacities(n_individual_compositions_and_phases, numbers::signaling_nan<double>()),
      compressibilities(n_individual_compositions_and_phases, numbers::signaling_nan<double>()),
      entropy_derivative_pressure(n_individual_compositions_and_phases, numbers::signaling_nan<double>()),
      entropy_derivative_temperature(n_individual_compositions_and_phases, numbers::signaling_nan<double>())
    {}



    template <int dim>
    void
    reaction_progress_modify_equation_of_state_outputs(const std::vector<double> &reaction_progress_values,
                                                       const std::vector<unsigned int> &reaction_progress_mapping,
                                                       const std::vector<unsigned int> &n_phase_transitions_per_composition,
                                                       EquationOfStateOutputs<dim> &eos_outputs_all_phases)
    {
      for (unsigned int c=0; c<n_phase_transitions_per_composition.size(); ++c)
        {
          MaterialModel::MaterialUtilities::reaction_progress_modify_values(reaction_progress_values, reaction_progress_mapping, n_phase_transitions_per_composition, eos_outputs_all_phases.densities, c);
          MaterialModel::MaterialUtilities::reaction_progress_modify_values(reaction_progress_values, reaction_progress_mapping, n_phase_transitions_per_composition, eos_outputs_all_phases.thermal_expansion_coefficients, c);
          MaterialModel::MaterialUtilities::reaction_progress_modify_values(reaction_progress_values, reaction_progress_mapping, n_phase_transitions_per_composition, eos_outputs_all_phases.specific_heat_capacities, c);
          MaterialModel::MaterialUtilities::reaction_progress_modify_values(reaction_progress_values, reaction_progress_mapping, n_phase_transitions_per_composition, eos_outputs_all_phases.compressibilities, c);
          MaterialModel::MaterialUtilities::reaction_progress_modify_values(reaction_progress_values, reaction_progress_mapping, n_phase_transitions_per_composition, eos_outputs_all_phases.entropy_derivative_pressure, c);
          MaterialModel::MaterialUtilities::reaction_progress_modify_values(reaction_progress_values, reaction_progress_mapping, n_phase_transitions_per_composition, eos_outputs_all_phases.entropy_derivative_temperature, c);
        }
    }



    template <int dim>
    void
    phase_average_equation_of_state_outputs(const EquationOfStateOutputs<dim> &eos_outputs_all_phases,
                                            const std::vector<double> &phase_function_values,
                                            const std::vector<unsigned int> &n_phase_transitions_per_composition,
                                            EquationOfStateOutputs<dim> &eos_outputs)
    {
      for (unsigned int c=0; c<eos_outputs.densities.size(); ++c)
        {
          eos_outputs.densities[c] =
            MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition, eos_outputs_all_phases.densities, c);
          eos_outputs.thermal_expansion_coefficients[c] =
            MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition, eos_outputs_all_phases.thermal_expansion_coefficients, c);
          eos_outputs.specific_heat_capacities[c] =
            MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition, eos_outputs_all_phases.specific_heat_capacities, c);
          eos_outputs.compressibilities[c] =
            MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition, eos_outputs_all_phases.compressibilities, c);
          eos_outputs.entropy_derivative_pressure[c] =
            MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition, eos_outputs_all_phases.entropy_derivative_pressure, c);
          eos_outputs.entropy_derivative_temperature[c] =
            MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition, eos_outputs_all_phases.entropy_derivative_temperature, c);
        }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  template struct EquationOfStateOutputs<dim>; \
  template void reaction_progress_modify_equation_of_state_outputs<dim>(const std::vector<double> &reaction_progress_values, \
                                                                        const std::vector<unsigned int> &reaction_progress_mapping, \
                                                                        const std::vector<unsigned int> &n_phase_transitions_per_composition, \
                                                                        EquationOfStateOutputs<dim> &); \
  template void phase_average_equation_of_state_outputs<dim> (const EquationOfStateOutputs<dim> &, \
                                                              const std::vector<double> &phase_function_values, \
                                                              const std::vector<unsigned int> &n_phase_transitions_per_composition, \
                                                              EquationOfStateOutputs<dim> &);

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
