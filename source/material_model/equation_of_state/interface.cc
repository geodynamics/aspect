/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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
  template void phase_average_equation_of_state_outputs<dim> (const EquationOfStateOutputs<dim> &, \
                                                              const std::vector<double> &phase_function_values, \
                                                              const std::vector<unsigned int> &n_phase_transitions_per_composition, \
                                                              EquationOfStateOutputs<dim> &);

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
