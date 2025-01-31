/*
  Copyright (C) 2022 - 2023 by the authors of the ASPECT code.

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
#include <aspect/global.h>
#include <aspect/simulator_signals.h>
#include <aspect/boundary_temperature/interface.h>

namespace aspect
{
  // Global variables (to be set by parameters)
  unsigned int switch_step;
  bool switched;

  /**
   * Declare additional parameters.
   */
  void declare_parameters(const unsigned int,
                          ParameterHandler &prm)
  {
    prm.declare_entry("Switch step", "0",
                      Patterns::Integer(0),
                      "Switch CFL at the timestep given.");
  }

  template <int dim>
  void parse_parameters(const Parameters<dim>,
                        ParameterHandler &prm)
  {
    switch_step = prm.get_integer("Switch step");
    switched = false;
  }

  template <int dim>
  void on_start_timestep (const SimulatorAccess<dim> &simulator_access)
  {
    simulator_access.get_pcout() << "Signal start_timestep triggered!" << std::endl;
    if (simulator_access.get_timestep_number() != numbers::invalid_unsigned_int
        &&
        simulator_access.get_timestep_number() >= switch_step
        &&
        !switched )
      {
        simulator_access.get_pcout() << "Reducing CFL number!" << std::endl;
        const_cast<Parameters<dim>&>(simulator_access.get_parameters()).CFL_number *= 0.5;

        switched = true;
      }
  }


  template <int dim>
  void on_edit_parameters (const SimulatorAccess<dim> &simulator_access,
                           Parameters<dim> &/*parameters*/)
  {
    simulator_access.get_pcout()<<"Signal edit_parameters triggered!"<<std::endl;
  }

  // Connect declare_parameters and parse_parameters to appropriate signals.
  void parameter_connector ()
  {
    SimulatorSignals<2>::declare_additional_parameters.connect (&declare_parameters);
    SimulatorSignals<3>::declare_additional_parameters.connect (&declare_parameters);

    SimulatorSignals<2>::parse_additional_parameters.connect (&parse_parameters<2>);
    SimulatorSignals<3>::parse_additional_parameters.connect (&parse_parameters<3>);
  }

  template <int dim>
  void signal_connector (SimulatorSignals<dim> &signals)
  {
    signals.edit_parameters_pre_setup_dofs.connect (&on_edit_parameters<dim>);
    signals.start_timestep.connect(&on_start_timestep<dim>);
  }

  // Tell ASPECT to send signals to the connector functions
  ASPECT_REGISTER_SIGNALS_PARAMETER_CONNECTOR(parameter_connector)
  ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>, signal_connector<3>)
}
