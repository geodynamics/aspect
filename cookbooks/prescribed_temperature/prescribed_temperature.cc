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

namespace aspect
{
  using namespace dealii;

  bool prescribe_internal_temperature;
  /**
   * Declare additional parameters.
   */
//  void declare_parameters(const unsigned int dim,
//                          ParameterHandler &prm)
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
  }



  template <int dim>
  void parse_parameters(const Parameters<dim> &,
                        ParameterHandler &prm)
  {
    prescribe_internal_temperature = prm.get_bool ("Prescribe internal temperature");
  }



  /**
   * This function is called by a signal which is triggered after the other constraints
   * have been calculated. This enables us to define additional constraints in the mass
   * matrix on any arbitrary degree of freedom in the model space.
   */
  template <int dim>
  void constrain_internal_temperature (const SimulatorAccess<dim> &simulator_access,
                                       AffineConstraints<double> &current_constraints)
  {
    if (prescribe_internal_temperature)
      {
        const std::vector< Point<dim>> points = aspect::Utilities::get_unit_support_points(simulator_access);
        const Quadrature<dim> quadrature (points);
        FEValues<dim> fe_values (simulator_access.get_fe(), quadrature, update_quadrature_points);
        typename DoFHandler<dim>::active_cell_iterator cell;
        MaterialModel::MaterialModelInputs<dim> in(quadrature.size(), simulator_access.introspection().n_compositional_fields);
        double num_constraints = 0;
        // Loop over all cells
        //       const unsigned int fixed_idx = simulator_access.introspection().compositional_index_for_name("Lithosphere");

        for (cell = simulator_access.get_dof_handler().begin_active();
             cell != simulator_access.get_dof_handler().end();
             ++cell)

          if (! cell->is_artificial())
            {
              fe_values.reinit (cell);
              // in.reint(cell, +args) in.reinit(fe_values, cell, this->introspection(), this->get_solution());
              std::vector<types::global_dof_index> local_dof_indices(simulator_access.get_fe().dofs_per_cell);
              cell->get_dof_indices (local_dof_indices);
              for (unsigned int q=0; q<quadrature.size(); q++)
                {
                  // If it's okay to constrain this DOF
                  if (current_constraints.can_store_line(local_dof_indices[q]) &&
                      !current_constraints.is_constrained(local_dof_indices[q]))
                    {

                    //  if (in.composition[fixed_idx][q] >= 0.5)
                      // Get the temperature component index
                      const unsigned int c_idx = simulator_access.get_fe().system_to_component_index(q).first;
                      // If we're on one of the temperature DOFs
                      if ((c_idx == simulator_access.introspection().component_indices.temperature))
                        {
                          const double fixed_temp = 2500;
                          current_constraints.add_line (local_dof_indices[q]);
                          current_constraints.set_inhomogeneity (local_dof_indices[q], fixed_temp);
                          num_constraints += 1;
                        }
                    }
                }
            }
      std::cout << num_constraints << std::endl;
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
