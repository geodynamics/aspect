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
  void post_nonlinear_solver (const SolverControl &nonlinear_solver_control)
  {
    const bool success = nonlinear_solver_control.last_check() == SolverControl::success;
    std::cout << "\nnumber of nonlinear iterations: " << nonlinear_solver_control.last_step() << ". State: " << success << ".\n";
  }


  template <int dim>
  void signal_connector (SimulatorSignals<dim> &signals)
  {
    std::cout << "Connecting signals" << std::endl;
    signals.post_nonlinear_solver.connect (&post_nonlinear_solver<dim>);
  }


  ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                    signal_connector<3>)
}
