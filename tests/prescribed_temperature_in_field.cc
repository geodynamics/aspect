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

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/fe/fe_values.h>
#include <aspect/global.h>
#include <aspect/utilities.h>
#include <aspect/simulator_signals.h>
#include <aspect/initial_temperature/interface.h>

namespace aspect
{
  // Declare and parse additional parameters
  bool prescribe_internal_temperature;
  std::vector <std::string> fixed_compositional_fields;
  double max_isotherm;


  void declare_parameters(const unsigned int /*dim*/,
                          ParameterHandler &prm)
  {
    prm.enter_subsection ("Prescribed internal temperature model");
    {
      prm.declare_entry ("Prescribe internal temperature", "false",
                         Patterns::Bool (),
                         "Whether or not to use any prescribed internal temperatures. "
                         "Locations in which to prescribe temperature are defined "
                         "based on the values of the compositional fields specified by "
                         "the parameter ``Names of compositional fields with fixed temperature'' "
                         "and the temperature is fixed to its initial state. Indicators are evaluated "
                         "separately at each support point."
                        );
      prm.declare_entry ("Maximum fixed temperature isosurface", "0",
                         Patterns::Double (),
                         "The maximum temperature that will remain fixed. Temperatures above this "
                         "value will be treated normally."
                        );
      prm.declare_entry ("Names of compositional fields with fixed temperature", "",
                         Patterns::List(Patterns::Anything ()),
                         "The names of compositional fields that determine the location where the "
                         "temperature is fixed to the values it was initialized with. The temperature "
                         "will be fixed if any of the fields listed here is above 0.5."
                        );
    }
    prm.leave_subsection ();
  }



  template <int dim>
  void parse_parameters(const Parameters<dim> &,
                        ParameterHandler &prm)
  {
    prm.enter_subsection ("Prescribed internal temperature model");
    {
      prescribe_internal_temperature = prm.get_bool ("Prescribe internal temperature");
      fixed_compositional_fields = Utilities::split_string_list (prm.get("Names of compositional fields with fixed temperature"));
      max_isotherm = prm.get_double ("Maximum fixed temperature isosurface");
    }
    prm.leave_subsection ();
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
    // Save access to the initial temperature manager the first time
    // we get here so that we can access it past the first time step
    // as well.
    static std::shared_ptr<const aspect::InitialTemperature::Manager<dim>> initial_temperature_manager;
    if (initial_temperature_manager == nullptr)
      initial_temperature_manager = simulator_access.get_initial_temperature_manager_pointer();

    if (prescribe_internal_temperature)
      {
        const std::vector<Point<dim>> points = aspect::Utilities::get_unit_support_points(simulator_access);
        const Quadrature<dim> quadrature (points);
        FEValues<dim> fe_values (simulator_access.get_fe(), quadrature, update_quadrature_points | update_values | update_gradients);
        typename DoFHandler<dim>::active_cell_iterator cell;
        MaterialModel::MaterialModelInputs<dim> in(quadrature.size(), simulator_access.introspection().n_compositional_fields);

        // Loop over all cells
        for (cell = simulator_access.get_dof_handler().begin_active();
             cell != simulator_access.get_dof_handler().end();
             ++cell)

          if (! cell->is_artificial())
            {
              fe_values.reinit (cell);
              in.reinit(fe_values, cell, simulator_access.introspection(), simulator_access.get_solution());

              std::vector<types::global_dof_index> local_dof_indices(simulator_access.get_fe().dofs_per_cell);
              cell->get_dof_indices (local_dof_indices);
              for (unsigned int q=0; q<quadrature.size(); q++)
                {

                  // If it's okay to constrain this DOF
                  if (current_constraints.can_store_line(local_dof_indices[q]) &&
                      !current_constraints.is_constrained(local_dof_indices[q]))
                    {

                      bool constrain_temperature = false;
                      for (unsigned int c=0; c < fixed_compositional_fields.size(); ++c)
                        {

                          const unsigned int composition_idx = simulator_access.introspection().compositional_index_for_name(fixed_compositional_fields[c]);
                          if (in.composition[q][composition_idx] >= 0.5 || in.temperature[q] <= max_isotherm)
                            {
                              constrain_temperature = true;
                            }
                        }

                      if (constrain_temperature)
                        {
                          // Get the temperature component index
                          const unsigned int c_idx = simulator_access.get_fe().system_to_component_index(q).first;

                          // If we're on one of the temperature DOFs
                          if (c_idx == simulator_access.introspection().component_indices.temperature)
                            {

                              // Set the temperature to be equal to the initial temperature
                              in.temperature[q] = initial_temperature_manager->initial_temperature(in.position[q]);
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
