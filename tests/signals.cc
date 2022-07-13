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

// make sure we can include deal.II and aspect files
#include <aspect/simulator_signals.h>

#include <iostream>

using namespace aspect;


// create a function that is run upon loading the plugin
// when declaring parameters, and that produces some output
void declare_parameters(const unsigned int dim,
                        ParameterHandler &prm)
{
  std::cout << "declaring parameters" << std::endl;
  prm.declare_entry("abc", "42", Patterns::Integer(21,43));
}


// same for parsing parameters
template <int dim>
void parse_parameters(const Parameters<dim> &parameters,
                      ParameterHandler &prm)
{
  std::cout << "parsing parameters: abc="
            << prm.get("abc") << std::endl;
  Assert (prm.get_integer ("abc") == 21, ExcInternalError());
}


void parameter_connector ()
{
  SimulatorSignals<2>::declare_additional_parameters.connect (&declare_parameters);
  SimulatorSignals<3>::declare_additional_parameters.connect (&declare_parameters);

  SimulatorSignals<2>::parse_additional_parameters.connect (&parse_parameters<2>);
  SimulatorSignals<3>::parse_additional_parameters.connect (&parse_parameters<3>);
}


ASPECT_REGISTER_SIGNALS_PARAMETER_CONNECTOR(parameter_connector)
