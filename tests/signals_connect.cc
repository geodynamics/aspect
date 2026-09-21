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

using namespace aspect;


template <int dim>
void post_constraints_creation (const SimulatorAccess<dim> &,
                                AffineConstraints<double> &)
{
  std::cout << std::endl
            << "*** constraints signal called ***"
            << std::endl
            << std::endl;
}


template <int dim>
void signal_connector (SimulatorSignals<dim> &signals)
{
  std::cout << "Connecting signals" << std::endl;
  signals.post_constraints_creation.connect (&post_constraints_creation<dim>);
}


ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                  signal_connector<3>)
