/*
  Copyright (C) 2024 - 2024 by the authors of the ASPECT code.

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

#include "common.h"

#include <aspect/parameters.h>
#include <aspect/simulator.h>
#include <aspect/introspection.h>

TEST_CASE("Introspection::basic")
{
  using namespace aspect;
  dealii::ParameterHandler prm;
  Simulator<2>::declare_parameters(prm);


  prm.set("Output directory", "");

  prm.enter_subsection("Compositional fields");
  prm.set("Number of fields", "3");
  prm.set("Names of fields", "a,b,c");
  prm.leave_subsection();

  Parameters<2> parameters(prm, MPI_COMM_WORLD);

  const std::vector<VariableDeclaration<2>> variables = construct_default_variables (parameters);
  Introspection<2> introspection (variables, parameters);

  dealii::FESystem<2> fe(introspection.get_fes(), introspection.get_multiplicities());
  CHECK(fe.get_name() == "FESystem<2>[FE_Q<2>(2)^2-FE_Q<2>(1)-FE_Q<2>(2)-FE_Q<2>(2)^3]");

  CHECK(introspection.n_components == 7);
  CHECK(introspection.n_compositional_fields == 3);
  CHECK(introspection.n_blocks == 6);

  CHECK(introspection.component_indices.velocities[0] == 0);
  CHECK(introspection.component_indices.velocities[1] == 1);
  CHECK(introspection.component_indices.pressure == 2);
  CHECK(introspection.component_indices.temperature == 3);
  CHECK(introspection.component_indices.compositional_fields[0] == 4);

  CHECK(introspection.block_indices.velocities == 0);
  CHECK(introspection.block_indices.pressure == 1);
  CHECK(introspection.block_indices.temperature == 2);
  CHECK(introspection.block_indices.compositional_fields.size() == 3);
  CHECK(introspection.block_indices.compositional_fields[0] == 3);
  CHECK(introspection.block_indices.compositional_fields[1] == 4);
  CHECK(introspection.block_indices.compositional_fields[2] == 5);

  // sparsity pattern block index:
  CHECK(introspection.block_indices.compositional_field_sparsity_pattern.size() == 3);
  CHECK(introspection.block_indices.compositional_field_sparsity_pattern[0] == 3);
  CHECK(introspection.block_indices.compositional_field_sparsity_pattern[1] == 3);
  CHECK(introspection.block_indices.compositional_field_sparsity_pattern[2] == 3);

  CHECK(introspection.base_elements.velocities == 0);
  CHECK(introspection.base_elements.pressure == 1);
  CHECK(introspection.base_elements.temperature == 2);
  // All compositional fields have the same FE:
  CHECK(introspection.base_elements.compositional_fields == std::vector<unsigned int>({3,3,3}));
  CHECK(introspection.get_composition_base_element_indices() == std::vector<unsigned int>({3}));
  CHECK(introspection.get_compositional_field_indices_with_base_element(3) == std::vector<unsigned int>({0,1,2}));
}

TEST_CASE("Introspection::different-composition-types")
{
  using namespace aspect;
  dealii::ParameterHandler prm;
  Simulator<2>::declare_parameters(prm);


  prm.set("Output directory", "");

  prm.enter_subsection("Compositional fields");
  prm.set("Number of fields", "3");
  prm.set("Names of fields", "a,b,c");
  prm.leave_subsection();
  prm.enter_subsection("Discretization");
  prm.set("Use discontinuous temperature discretization", "false");
  prm.set("Use discontinuous composition discretization", "true,false,false");
  prm.leave_subsection();

  Parameters<2> parameters(prm, MPI_COMM_WORLD);

  const std::vector<VariableDeclaration<2>> variables = construct_default_variables (parameters);
  Introspection<2> introspection (variables, parameters);

  dealii::FESystem<2> fe(introspection.get_fes(), introspection.get_multiplicities());
  CHECK(fe.get_name() == "FESystem<2>[FE_Q<2>(2)^2-FE_Q<2>(1)-FE_Q<2>(2)-FE_DGQ<2>(2)-FE_Q<2>(2)^2]");

  CHECK(introspection.n_components == 7);
  CHECK(introspection.n_compositional_fields == 3);
  CHECK(introspection.n_blocks == 6);

  CHECK(introspection.component_indices.velocities[0] == 0);
  CHECK(introspection.component_indices.velocities[1] == 1);
  CHECK(introspection.component_indices.pressure == 2);
  CHECK(introspection.component_indices.temperature == 3);
  CHECK(introspection.component_indices.compositional_fields[0] == 4);

  CHECK(introspection.block_indices.velocities == 0);
  CHECK(introspection.block_indices.pressure == 1);
  CHECK(introspection.block_indices.temperature == 2);
  CHECK(introspection.block_indices.compositional_fields[0] == 3);
  CHECK(introspection.block_indices.compositional_fields[1] == 4);
  CHECK(introspection.block_indices.compositional_fields[2] == 5);

  // sparsity pattern block index:
  CHECK(introspection.block_indices.compositional_field_sparsity_pattern[0] == 3);
  CHECK(introspection.block_indices.compositional_field_sparsity_pattern[1] == 4);
  CHECK(introspection.block_indices.compositional_field_sparsity_pattern[2] == 4);

  // base elements:
  CHECK(introspection.base_elements.compositional_fields == std::vector<unsigned int>({3,4,4}));
  CHECK(introspection.get_composition_base_element_indices() == std::vector<unsigned int>({3,4}));
  CHECK(introspection.get_compositional_field_indices_with_base_element(3) == std::vector<unsigned int>({0}));
  CHECK(introspection.get_compositional_field_indices_with_base_element(4) == std::vector<unsigned int>({1,2}));
}
