// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2018-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include "common.h"
#include <aspect/utilities.h>

// A test that simply verifies that the additional material
// model outputs can be accessed.

#include <aspect/simulator.h>
#include <deal.II/grid/tria.h>
#include <aspect/material_model/simple.h>
#include <aspect/simulator_access.h>

#include <iostream>

namespace
{
  using namespace dealii;
  using namespace aspect;

  template <int dim>
  class AdditionalOutputs1 : public MaterialModel::AdditionalMaterialOutputs<dim>
  {
    public:
      AdditionalOutputs1 (const unsigned int n_points,
                          const unsigned int /*n_comp*/)
      {
        additional_material_output1.resize(n_points);
      }

      std::vector<double> additional_material_output1;
  };


  template <int dim>
  class Material1 : public MaterialModel::Simple<dim>
  {
    public:

      void evaluate(const MaterialModel::MaterialModelInputs<dim> &/*in*/,
                    MaterialModel::MaterialModelOutputs<dim> &out) const override
      {
        const auto additional = out.template get_additional_output_object<AdditionalOutputs1<dim>>();
        additional->additional_material_output1[0] = 42.0;
      }
  };
}


TEST_CASE("AdditionalOutputs works")
{
  const int dim=2;

  using namespace aspect::MaterialModel;
  MaterialModelInputs<dim> in(1,1);
  MaterialModelOutputs<dim> out(1,1);
  in.requested_properties = MaterialProperties::additional_outputs;


  REQUIRE(out.has_additional_output_object<AdditionalOutputs1<dim>>() == false);

  out.additional_outputs.push_back(std::make_unique<AdditionalOutputs1<dim>> (1, 1));

  REQUIRE(out.has_additional_output_object<AdditionalOutputs1<dim>>());

  Material1<dim> mat;
  mat.evaluate(in, out);

  REQUIRE(out.get_additional_output_object<AdditionalOutputs1<dim>>()->additional_material_output1[0] == 42.0);

  // test const version of get_additional_output:
  {
    const MaterialModelOutputs<dim> &const_out = out;
    REQUIRE(const_out.has_additional_output_object<AdditionalOutputs1<dim>>());
    const auto a = const_out.get_additional_output_object<AdditionalOutputs1<dim>>();
    REQUIRE(a != nullptr);
  }
}
