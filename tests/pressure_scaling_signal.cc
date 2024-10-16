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
