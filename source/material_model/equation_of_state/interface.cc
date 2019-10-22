/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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
    EquationOfStateOutputs<dim>::EquationOfStateOutputs(const unsigned int n_compositions)
      :
      densities(n_compositions, numbers::signaling_nan<double>()),
      thermal_expansion_coefficients(n_compositions, numbers::signaling_nan<double>()),
      specific_heat_capacities(n_compositions, numbers::signaling_nan<double>()),
      compressibilities(n_compositions, numbers::signaling_nan<double>()),
      entropy_derivative_pressure(n_compositions, numbers::signaling_nan<double>()),
      entropy_derivative_temperature(n_compositions, numbers::signaling_nan<double>())
    {}

    template <int dim>
    void
    compute_equation_of_state_phase_transitions(const EquationOfStateOutputs<dim> &eos_outputs_all_phases,
                                                const MaterialUtilities::PhaseFunction<dim> &phase_function,
                                                MaterialUtilities::PhaseFunctionInputs<dim> &phase_in,
                                                EquationOfStateOutputs<dim> &eos_outputs)
    {
      unsigned int j=0;
      for (unsigned int c=0; c<eos_outputs.densities.size(); ++c)
        {
          Assert(j<eos_outputs_all_phases.densities.size(),
                 ExcInternalError());

          eos_outputs.densities[c] = eos_outputs_all_phases.densities[j];
          eos_outputs.thermal_expansion_coefficients[c] = eos_outputs_all_phases.thermal_expansion_coefficients[j];
          eos_outputs.specific_heat_capacities[c] = eos_outputs_all_phases.specific_heat_capacities[j];
          eos_outputs.compressibilities[c] = eos_outputs_all_phases.compressibilities[j];
          eos_outputs.entropy_derivative_pressure[c] = eos_outputs_all_phases.entropy_derivative_pressure[j];
          eos_outputs.entropy_derivative_temperature[c] = eos_outputs_all_phases.entropy_derivative_temperature[j];

          const unsigned int base_index = j;
          ++j;

          for (unsigned int p=0; p<phase_function.n_phase_transitions_for_composition(c); ++p)
            {
              Assert(j<eos_outputs_all_phases.densities.size(),
                     ExcInternalError());

              phase_in.phase_index = j-c-1;
              eos_outputs.densities[c] += phase_function.compute_value(phase_in) * (eos_outputs_all_phases.densities[j]-eos_outputs_all_phases.densities[base_index]);
              eos_outputs.thermal_expansion_coefficients[c] += phase_function.compute_value(phase_in) * (eos_outputs_all_phases.thermal_expansion_coefficients[j]-eos_outputs_all_phases.thermal_expansion_coefficients[base_index]);
              eos_outputs.specific_heat_capacities[c] += phase_function.compute_value(phase_in) * (eos_outputs_all_phases.specific_heat_capacities[j]-eos_outputs_all_phases.specific_heat_capacities[base_index]);

              ++j;
            }
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
  template void compute_equation_of_state_phase_transitions<dim> (const EquationOfStateOutputs<dim> &, \
                                                                  const MaterialUtilities::PhaseFunction<dim> &, \
                                                                  MaterialUtilities::PhaseFunctionInputs<dim> &, \
                                                                  EquationOfStateOutputs<dim> &);

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}
