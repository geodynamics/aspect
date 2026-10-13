// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2022-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

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

    if (current_constraints.can_store_line(first_pressure_dof))
      {
        current_constraints.add_line(first_pressure_dof);
        current_constraints.set_inhomogeneity (first_pressure_dof, value);

        std::cout << "adding constraint idx= " << first_pressure_dof
                  << " to value= " << value
                  << std::endl;
      }
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
