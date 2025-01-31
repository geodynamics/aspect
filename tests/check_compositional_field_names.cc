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
#include <aspect/parameters.h>

template <int dim>
void f(const aspect::SimulatorAccess<dim> &simulator_access,
       aspect::Assemblers::Manager<dim> &)
{
  // This function tests whether the compositional field types are correctly
  // processed.

  using namespace aspect::MaterialModel;

  aspect::ParameterHandler prm;

  const std::vector<std::string> c_names = simulator_access.introspection().get_composition_names();
  const std::vector<typename aspect::Parameters<dim>::CompositionalFieldDescription> descriptions = simulator_access.introspection().get_composition_descriptions();

  for (unsigned int i=0; i<simulator_access.introspection().n_compositional_fields; ++i)
    {
      if (descriptions[i].type == aspect::Parameters<dim>::CompositionalFieldDescription::chemical_composition)
        std::cout << c_names[i] << " is of type chemical composition" << std::endl;
      if (descriptions[i].type == aspect::Parameters<dim>::CompositionalFieldDescription::grain_size)
        std::cout << c_names[i] << " is of type grain size" << std::endl;
      if (descriptions[i].type == aspect::Parameters<dim>::CompositionalFieldDescription::porosity)
        std::cout << c_names[i] << " is of type porosity" << std::endl;
      if (descriptions[i].type == aspect::Parameters<dim>::CompositionalFieldDescription::generic)
        std::cout << c_names[i] << " is of type generic" << std::endl;
      if (descriptions[i].type == aspect::Parameters<dim>::CompositionalFieldDescription::stress)
        std::cout << c_names[i] << " is of type stress" << std::endl;
      if (descriptions[i].type == aspect::Parameters<dim>::CompositionalFieldDescription::unspecified)
        std::cout << c_names[i] << " is of type unspecified" << std::endl;
    }

  exit(0);

}

template <>
void f(const aspect::SimulatorAccess<2> &,
       aspect::Assemblers::Manager<2> &)
{
  AssertThrow(false,dealii::ExcInternalError());
}

template <int dim>
void signal_connector (aspect::SimulatorSignals<dim> &signals)
{
  std::cout << "* Connecting signals" << std::endl;
  signals.set_assemblers.connect (std::bind(&f<dim>,
                                            std::placeholders::_1,
                                            std::placeholders::_2));
}

ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                  signal_connector<3>)
