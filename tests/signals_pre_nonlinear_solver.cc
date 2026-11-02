// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2022-2026 by the ASPECT authors.
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

using namespace aspect;


template <int dim>
void pre_nonlinear_solver (const SimulatorAccess<dim> &simulator_access)
{
  simulator_access.get_pcout() << "Signal was triggered" << std::endl;
}


template <int dim>
void signal_connector (SimulatorSignals<dim> &signals)
{
  signals.pre_nonlinear_solver.connect (&pre_nonlinear_solver<dim>);
}


ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                  signal_connector<3>)
