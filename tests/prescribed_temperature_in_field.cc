/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/fe/fe_values.h>
#include <aspect/global.h>
#include <aspect/utilities.h>
#include <aspect/simulator_signals.h>
#include <aspect/initial_temperature/interface.h>

namespace aspect
{
  using namespace dealii;

  // Declare and parse additional parameters
  bool prescribe_internal_temperature;
  unsigned int fixed_composition_idx;
  void declare_parameters(const unsigned int dim,
                          ParameterHandler &prm)
  {
    prm.declare_entry ("Prescribe internal temperature", "false",
                       Patterns::Bool (),
                       "Whether or not to use any prescribed internal temperature. "
                       "Locations in which to prescribe temperature are defined "
                       "in section ``Prescribed temperature/Indicator function'' "
                       "and the temperature are defined in section ``Prescribed "
                       "temperature/Temperature function''. Indicators are evaluated "
                       "at the center of each cell, and all DOFs associated with "
                       "the specified temperature at the indicated cells "
                       "are constrained."
                      );
    prm.declare_entry ("Index of fixed compositional field", "0",
                       Patterns::Integer (0),
                       "The index of the compositional field that will have its "
                       "temperature fixed to the values it was initialized with."
                      );
  }

  template <int dim>
  void parse_parameters(const Parameters<dim> &,
                        ParameterHandler &prm)
  {
    prescribe_internal_temperature = prm.get_bool ("Prescribe internal temperature");
    fixed_composition_idx = prm.get_integer ("Index of fixed compositional field");
  }

  /**
   * This function is called by a signal which is triggered after the other constraints
   * have been calculated. This enables the definition of additional constraints in the mass
   * matrix on any arbitrary degree of freedom in the model space.
   */
  template <int dim>
  void constrain_internal_temperature (const SimulatorAccess<dim> &simulator_access,
                                       AffineConstraints<double> &current_constraints)
  {
    if (prescribe_internal_temperature)
      {
        const std::vector<Point<dim>> points = aspect::Utilities::get_unit_support_points(simulator_access);
        const Quadrature<dim> quadrature (points);
        FEValues<dim> fe_values (simulator_access.get_fe(), quadrature, update_quadrature_points | update_values | update_gradients); // update_values
        typename DoFHandler<dim>::active_cell_iterator cell;
        MaterialModel::MaterialModelInputs<dim> in(quadrature.size(), simulator_access.introspection().n_compositional_fields);

        // Loop over all cells
        for (cell = simulator_access.get_dof_handler().begin_active();
            cell != simulator_access.get_dof_handler().end();
            ++cell)

          if (! cell->is_artificial())
            {
              fe_values.reinit (cell);
              in.reinit(fe_values, cell, simulator_access.introspection(), simulator_access.get_solution(), false);
              
              std::vector<types::global_dof_index> local_dof_indices(simulator_access.get_fe().dofs_per_cell);
              cell->get_dof_indices (local_dof_indices);
              for (unsigned int q=0; q<quadrature.size(); q++)
                {
                  // If it's okay to constrain this DOF
                  if (current_constraints.can_store_line(local_dof_indices[q]) &&
                      !current_constraints.is_constrained(local_dof_indices[q]))
                    {
                      // Make sure we're in the desired compositional field
                      if (in.composition[q][fixed_composition_idx] >= 0.5)
                        {
                          // Get the temperature component index
                          const unsigned int c_idx = simulator_access.get_fe().system_to_component_index(q).first;
                          // If we're on one of the temperature DOFs
                          if (c_idx == simulator_access.introspection().component_indices.temperature)
                            {
                              // Set the temperature to be equal to the initial temperature
                              in.temperature[q] = simulator_access.get_initial_temperature_manager().initial_temperature(in.position[q]);
                              // Update the constraints
                              current_constraints.add_line (local_dof_indices[q]);
                              current_constraints.set_inhomogeneity (local_dof_indices[q], in.temperature[q]);
                            }
                        }
                    }
                }
            }
      }
  }

  // Connect declare_parameters and parse_parameters to appropriate signals.
  void parameter_connector ()
  {
    SimulatorSignals<2>::declare_additional_parameters.connect (&declare_parameters);
    SimulatorSignals<3>::declare_additional_parameters.connect (&declare_parameters);

    SimulatorSignals<2>::parse_additional_parameters.connect (&parse_parameters<2>);
    SimulatorSignals<3>::parse_additional_parameters.connect (&parse_parameters<3>);
  }

  // Connect constraints function to correct signal.
  template <int dim>
  void signal_connector (SimulatorSignals<dim> &signals)
  {
    signals.post_constraints_creation.connect (&constrain_internal_temperature<dim>);
  }

  // Tell ASPECT to send signals to the connector functions
  ASPECT_REGISTER_SIGNALS_PARAMETER_CONNECTOR(parameter_connector)
  ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>, signal_connector<3>)
}
