/*
  Copyright (C) 2022 - 2024 by the authors of the ASPECT code.

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
