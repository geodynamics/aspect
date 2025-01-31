/*
  Copyright (C) 2023 by the authors of the ASPECT code.

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
#include <aspect/material_model/utilities.h>

template <int dim>
void f(const aspect::SimulatorAccess<dim> &simulator_access,
       aspect::Assemblers::Manager<dim> &)
{
  // This function tests whether the compositional field types are correctly
  // processed.

  using namespace aspect::MaterialModel;

  aspect::ParameterHandler prm;

  // Fields 0 and 4 are chemical compositions
  // These are called Field 1 and Field 5
  const std::vector<unsigned int> indices = simulator_access.introspection().get_indices_for_fields_of_type(aspect::Parameters<dim>::CompositionalFieldDescription::chemical_composition);
  const std::vector<std::string> names = simulator_access.introspection().get_names_for_fields_of_type(aspect::Parameters<dim>::CompositionalFieldDescription::chemical_composition);

  const std::vector<double> field_values { 0.1, 0.3, 0.3, 0.3, 0.2, 0.3, 0.3};
  const std::vector<double> compositional_field_fractions = MaterialUtilities::compute_only_composition_fractions(field_values, indices);

  std::cout << "Chemical composition field #0 (the background field) has value " << compositional_field_fractions[0] << "." << std::endl;

  for (unsigned int i=0; i<3; ++i)
    std::cout << "Chemical composition field #" << i+1 << " is called " << names[i] << ", is in position "  << indices[i] << ", and has value " << compositional_field_fractions[i+1] << "." << std::endl;

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
