/*
  Copyright (C) 2018 - 2022 by the authors of the ASPECT code.

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
        AdditionalOutputs1<dim> *additional;

        additional = out.template get_additional_output<AdditionalOutputs1<dim>>();
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


  REQUIRE(out.get_additional_output<AdditionalOutputs1<dim>>() == NULL);

  out.additional_outputs.push_back(std::make_unique<AdditionalOutputs1<dim>> (1, 1));

  REQUIRE(out.get_additional_output<AdditionalOutputs1<dim>>() != NULL);

  Material1<dim> mat;
  mat.evaluate(in, out);

  REQUIRE(out.get_additional_output<AdditionalOutputs1<dim>>()->additional_material_output1[0] == 42.0);

  // test const version of get_additional_output:
  {
    const MaterialModelOutputs<dim> &const_out = out;
    REQUIRE(const_out.get_additional_output<AdditionalOutputs1<dim>>() != NULL);
    const AdditionalOutputs1<dim> *a = const_out.get_additional_output<AdditionalOutputs1<dim>>();
    REQUIRE(a != nullptr);
  }
}
