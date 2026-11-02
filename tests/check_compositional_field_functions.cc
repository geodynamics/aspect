// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2023-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/simulator_signals.h>
#include <aspect/introspection.h>
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
  const std::vector<unsigned int> indices = simulator_access.introspection().get_indices_for_fields_of_type(aspect::CompositionalFieldDescription::chemical_composition);
  const std::vector<std::string> names = simulator_access.introspection().get_names_for_fields_of_type(aspect::CompositionalFieldDescription::chemical_composition);

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
