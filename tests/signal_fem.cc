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

#include <aspect/simulator.h>
#include <deal.II/grid/tria.h>
#include <aspect/material_model/simple.h>
#include <aspect/simulator_access.h>
#include <aspect/simulator/assemblers/interface.h>

#include <deal.II/fe/fe_dgq.h>
#include <iostream>

namespace aspect
{

  template <int dim>
  void my_signal(std::vector<VariableDeclaration<dim>> &variables)
  {
    std::cout << "* signals.edit_finite_element_variables:" << std::endl;

    VariableDeclaration<dim> dummy("dummy",
                                   std::make_shared<FE_DGQ<dim>>(4),
                                   2,
                                   1);
    variables.insert(variables.begin()+3, dummy);

    for (unsigned int i=0; i<variables.size(); ++i)
      {
        std::cout << " name=" << variables[i].name
                  << " fe=" << variables[i].fe->get_name()
                  << " multiplicity=" << variables[i].multiplicity
                  << " n_blocks=" << variables[i].n_blocks
                  << " n_components=" << variables[i].n_components()
                  << std::endl;
      }

    std::cout << std::endl;
  }
}




template <int dim>
void signal_connector (aspect::SimulatorSignals<dim> &signals)
{
  std::cout << "* Connecting signals" << std::endl;
  signals.edit_finite_element_variables.connect(&aspect::my_signal<dim>);
}

ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                  signal_connector<3>)
