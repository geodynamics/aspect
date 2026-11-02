// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2022-2024 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include "../benchmarks/inclusion/inclusion.cc"

#include <aspect/simulator_signals.h>


double my_signal(const double pressure_scaling, const double reference_viscosity, const double length_scale)
{
  std::cout << "pressure_scaling = " << pressure_scaling
            << " reference_viscosity = " << reference_viscosity
            << " length_scale = " << length_scale << std::endl;
  // Now change it to something different:
  return 42.0;
}


template <int dim>
void signal_connector (aspect::SimulatorSignals<dim> &signals)
{
  std::cout << "* Connecting signals" << std::endl;
  signals.modify_pressure_scaling.connect(&my_signal);
}

ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                  signal_connector<3>)
