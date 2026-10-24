// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2022 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include <aspect/simulator_signals.h>
#include <aspect/simulator_access.h>

#include <iostream>

namespace aspect
{
  template <int dim>
  void post_stokes_solver (const SimulatorAccess<dim> &,
                           const unsigned int /*number_S_iterations*/,
                           const unsigned int /*number_A_iterations*/,
                           const SolverControl &/*solver_control_cheap*/,
                           const SolverControl &/*solver_control_expensive*/)
  {
    std::cout << "\npost_stokes_solver:\n";
  }


  template <int dim>
  void signal_connector (SimulatorSignals<dim> &signals)
  {
    std::cout << "Connecting signals" << std::endl;
    signals.post_stokes_solver.connect (&post_stokes_solver<dim>);
  }


  ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                    signal_connector<3>)
}
