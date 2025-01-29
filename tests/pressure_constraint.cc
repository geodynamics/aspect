/*
  Copyright (C) 2022 by the authors of the ASPECT code.

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
#include <aspect/simulator_signals.h>

namespace aspect
{
  template <int dim>
  void modify_constraints (const SimulatorAccess<dim> &simulator_access,
                           AffineConstraints<double> &current_constraints)
  {
    // Hack: the first pressure dof is only this easy to compute if we don't
    // use a direct solver or reorganize the blocks of the linear system in
    // some other way. Good enough for this test though!
    const types::global_dof_index first_pressure_dof =
      simulator_access.introspection().system_dofs_per_block[0];

    const double value = 100;

    current_constraints.add_line(first_pressure_dof);
    current_constraints.set_inhomogeneity (first_pressure_dof, value);

    std::cout << "adding constraint idx= " << first_pressure_dof
              << " to value= " << value
              << std::endl;
  }

  // Connect constraints function to correct signal.
  template <int dim>
  void signal_connector (SimulatorSignals<dim> &signals)
  {
    signals.post_constraints_creation.connect (&modify_constraints<dim>);
  }

  // Tell ASPECT to send signals to the connector functions
  ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>, signal_connector<3>)
}
