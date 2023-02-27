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

TEST_CASE("Utilities::weighted_p_norm_average")
{
  std::vector<double> weights = {1,1,2,2,3,3};
  std::vector<double> values = {6,5,4,3,2,1};
  std::vector<double> p_norms = {-1000,-2.5,-2,-1,0,1,2,2.5,3,4,1000};
  std::vector<double> expected = {1., 1.59237, 1.6974 , 1.98895, 2.38899, 2.83333, 3.24037, 3.41824, 3.57872, 3.85347, 6. };

  for (unsigned int i = 0; i < p_norms.size(); i++)
    {
      INFO("check i=" << i << ": ");
      REQUIRE(aspect::Utilities::weighted_p_norm_average(weights,values,p_norms[i]) == Approx(expected[i]));
    }

}

TEST_CASE("Utilities::AsciiDataLookup")
{
  using namespace dealii;

  //TODO: add support for setting data directly instead of relying on a file to load:
  aspect::Utilities::StructuredDataLookup<1> lookup(2 /*n_components*/, 1.0 /*scaling*/);
  lookup.load_file(ASPECT_SOURCE_DIR "/data/boundary-velocity/ascii-data/test/box_2d_left.0.txt", MPI_COMM_WORLD);

  INFO(lookup.get_data(Point<1>(330000./2.0),0));
  INFO(lookup.get_data(Point<1>(330000./2.0),1));
  INFO(lookup.get_gradients(Point<1>(330000./2.0),0));
  INFO(lookup.get_gradients(Point<1>(330000./2.0),1));

  REQUIRE(lookup.get_data(Point<1>(330000./2.0),0) == Approx(0.5));
  REQUIRE(lookup.get_data(Point<1>(330000./2.0),1) == Approx(0.0));
  REQUIRE(lookup.get_gradients(Point<1>(330000./2.0),0)[0] == Approx(-1.0/330000.));
  REQUIRE(lookup.get_gradients(Point<1>(330000./2.0),1)[0] == Approx(0.0));
}


TEST_CASE("Utilities::AsciiDataLookup manual dim=1")
{
  using namespace dealii;

  aspect::Utilities::StructuredDataLookup<1> lookup(2 /*n_components*/, 1.0 /*scaling*/);

  std::vector<std::string> column_names = {"a", "b"};
  Table<1,double> table(2);
  std::vector<Table<1,double>> raw_data(2, table);

  std::vector<std::vector<double>> coordinate_values(1, std::vector<double>({1.0, 2.0}));
  // c1:
  raw_data[0](0) = 0.0;
  raw_data[0](1) = 1.0;
  // c2:
  raw_data[1](0) = 5.0;
  raw_data[1](1) = 3.0;

  lookup.reinit(column_names, std::move(coordinate_values), std::move(raw_data),
                MPI_COMM_SELF, numbers::invalid_unsigned_int);

  INFO(lookup.get_data(Point<1>(1.5), 0));
  INFO(lookup.get_data(Point<1>(1.5), 1));

  REQUIRE(lookup.get_data(Point<1>(1.5),0) == Approx(0.5));
  REQUIRE(lookup.get_data(Point<1>(1.5),1) == Approx(4.0));

  REQUIRE(lookup.get_gradients(Point<1>(1.5),0)[0] == Approx(1.0));
  REQUIRE(lookup.get_gradients(Point<1>(1.5),1)[0] == Approx(-2.0));
}

TEST_CASE("Utilities::AsciiDataLookup manual dim=2")
{
  using namespace dealii;

  aspect::Utilities::StructuredDataLookup<2> lookup(1 /*n_components*/, 1.0 /*scaling*/);

  std::vector<std::string> column_names = {"topography"};
  std::vector<Table<2,double>> raw_data(1, Table<2,double>(3,3));
  std::vector<std::vector<double>> coordinate_values(2, std::vector<double>(3, 0.));

  // x:
  coordinate_values[0] = {0., 1., 3.};
  // y:
  coordinate_values[1] = {5., 6., 7.};
  // c1:
  raw_data[0](0,0) = 1.0;
  raw_data[0](1,0) = 2.0;
  raw_data[0](2,0) = 3.0;
  raw_data[0](0,1) = 4.0;
  raw_data[0](1,1) = 5.0;
  raw_data[0](2,1) = 6.0;
  raw_data[0](0,2) = 0.0;
  raw_data[0](1,2) = 0.0;
  raw_data[0](2,2) = 0.0;

  lookup.reinit(column_names, std::move(coordinate_values), std::move(raw_data),
                MPI_COMM_SELF, numbers::invalid_unsigned_int);

  REQUIRE(lookup.get_data(Point<2>(1.0,6.0),0) == Approx(5.0));
  REQUIRE(lookup.get_data(Point<2>(2.0,6.0),0) == Approx(5.5));
}

TEST_CASE("Utilities::AsciiDataLookup manual dim=2 equid")
{
  using namespace dealii;

  aspect::Utilities::StructuredDataLookup<2> lookup(1 /*n_components*/, 1.0 /*scaling*/);

  std::vector<std::string> column_names = {"topography"};
  std::vector<Table<2,double>> raw_data(1, Table<2,double>(3,3));
  std::vector<std::vector<double>> coordinate_values(2, std::vector<double>(3, 0.));

  // x:
  coordinate_values[0] = {0., 1., 2.};
  // y:
  coordinate_values[1] = {5., 6., 7.};
  // c1:
  raw_data[0](0,0) = 1.0;
  raw_data[0](1,0) = 2.0;
  raw_data[0](2,0) = 3.0;
  raw_data[0](0,1) = 4.0;
  raw_data[0](1,1) = 5.0;
  raw_data[0](2,1) = 6.0;
  raw_data[0](0,2) = 0.0;
  raw_data[0](1,2) = 0.0;
  raw_data[0](2,2) = 0.0;

  lookup.reinit(column_names, std::move(coordinate_values), std::move(raw_data),
                MPI_COMM_SELF, numbers::invalid_unsigned_int);

  REQUIRE(lookup.get_data(Point<2>(1.0,6.0),0) == Approx(5.0));
  REQUIRE(lookup.get_data(Point<2>(1.5,6.0),0) == Approx(5.5));
}

TEST_CASE("Random draw volume weighted average rotation matrix")
{
  std::vector<double> unsorted_volume_fractions = {2.,5.,1.,3.,6.,4.};
  std::vector<double> sorted_volume_fractions_ref = {1.,2.,3.,4.,5.,6.};
  const std::vector<std::size_t> permutation = aspect::Utilities::sort_permutation<double>(unsorted_volume_fractions);
  const std::vector<double> sorted_volume_fractions = aspect::Utilities::apply_permutation<double>(unsorted_volume_fractions,permutation);

  for (unsigned int i = 0; i < sorted_volume_fractions.size(); i++)
    {
      REQUIRE(sorted_volume_fractions[i] == Approx(sorted_volume_fractions_ref[i]));
    }

  const std::vector<dealii::Tensor<2,3>> unsorted_rotation_matrices =
  {
    dealii::Tensor<2,3>({{2.,2.,2.},{2.,2.,2.},{2.,2.,2.}}),
    dealii::Tensor<2,3>({{5.,5.,5.},{5.,5.,5.},{5.,5.,5.}}),
    dealii::Tensor<2,3>({{1.,1.,1.},{1.,1.,1.},{1.,1.,1.}}),
    dealii::Tensor<2,3>({{3.,3.,3.},{3.,3.,3.},{3.,3.,3.}}),
    dealii::Tensor<2,3>({{6.,6.,6.},{6.,6.,6.},{6.,6.,6.}}),
    dealii::Tensor<2,3>({{4.,4.,4.},{4.,4.,4.},{4.,4.,4.}})
  };

  const std::vector<dealii::Tensor<2,3>> sorted_rotation_matrices_ref =
  {
    dealii::Tensor<2,3>({{1.,1.,1.},{1.,1.,1.},{1.,1.,1.}}),
    dealii::Tensor<2,3>({{2.,2.,2.},{2.,2.,2.},{2.,2.,2.}}),
    dealii::Tensor<2,3>({{3.,3.,3.},{3.,3.,3.},{3.,3.,3.}}),
    dealii::Tensor<2,3>({{4.,4.,4.},{4.,4.,4.},{4.,4.,4.}}),
    dealii::Tensor<2,3>({{5.,5.,5.},{5.,5.,5.},{5.,5.,5.}}),
    dealii::Tensor<2,3>({{6.,6.,6.},{6.,6.,6.},{6.,6.,6.}})
  };
  const std::vector<dealii::Tensor<2,3>> sorted_rotation_matrices = aspect::Utilities::apply_permutation<dealii::Tensor<2,3>>(unsorted_rotation_matrices,permutation);
  for (unsigned int i = 0; i < sorted_rotation_matrices.size(); i++)
    {
      REQUIRE(sorted_rotation_matrices[i][0][0] == Approx(sorted_rotation_matrices[i][0][0]));
    }

  std::mt19937 random_number_generator;
  random_number_generator.seed(5);
  const std::vector<dealii::Tensor<2,3>> result = aspect::Utilities::random_draw_volume_weighting_rotation_matrices(unsorted_volume_fractions,
                                                   unsorted_rotation_matrices,
                                                   25,
                                                   random_number_generator);

  const std::vector<dealii::Tensor<2,3>> result_ref =
  {
    dealii::Tensor<2,3>({{2.,2.,2.},{2.,2.,2.},{2.,2.,2.}}),
    dealii::Tensor<2,3>({{6.,6.,6.},{6.,6.,6.},{6.,6.,6.}}),
    dealii::Tensor<2,3>({{4.,4.,4.},{4.,4.,4.},{4.,4.,4.}}),
    dealii::Tensor<2,3>({{6.,6.,6.},{6.,6.,6.},{6.,6.,6.}}),
    dealii::Tensor<2,3>({{2.,2.,2.},{2.,2.,2.},{2.,2.,2.}}),
    dealii::Tensor<2,3>({{4.,4.,4.},{4.,4.,4.},{4.,4.,4.}}),
    dealii::Tensor<2,3>({{4.,4.,4.},{4.,4.,4.},{4.,4.,4.}}),
    dealii::Tensor<2,3>({{5.,5.,5.},{5.,5.,5.},{5.,5.,5.}}),
    dealii::Tensor<2,3>({{6.,6.,6.},{6.,6.,6.},{6.,6.,6.}}),
    dealii::Tensor<2,3>({{6.,6.,6.},{6.,6.,6.},{6.,6.,6.}}),
    dealii::Tensor<2,3>({{5.,5.,5.},{5.,5.,5.},{5.,5.,5.}}),
    dealii::Tensor<2,3>({{6.,6.,6.},{6.,6.,6.},{6.,6.,6.}}),
    dealii::Tensor<2,3>({{3.,3.,3.},{3.,3.,3.},{3.,3.,3.}}),
    dealii::Tensor<2,3>({{6.,6.,6.},{6.,6.,6.},{6.,6.,6.}}),
    dealii::Tensor<2,3>({{2.,2.,2.},{2.,2.,2.},{2.,2.,2.}}),
    dealii::Tensor<2,3>({{3.,3.,3.},{3.,3.,3.},{3.,3.,3.}}),
    dealii::Tensor<2,3>({{2.,2.,2.},{2.,2.,2.},{2.,2.,2.}}),
    dealii::Tensor<2,3>({{6.,6.,6.},{6.,6.,6.},{6.,6.,6.}}),
    dealii::Tensor<2,3>({{5.,5.,5.},{5.,5.,5.},{5.,5.,5.}}),
    dealii::Tensor<2,3>({{5.,5.,5.},{5.,5.,5.},{5.,5.,5.}}),
    dealii::Tensor<2,3>({{6.,6.,6.},{6.,6.,6.},{6.,6.,6.}}),
    dealii::Tensor<2,3>({{1.,1.,1.},{1.,1.,1.},{1.,1.,1.}}),
    dealii::Tensor<2,3>({{2.,2.,2.},{2.,2.,2.},{2.,2.,2.}}),
    dealii::Tensor<2,3>({{6.,6.,6.},{6.,6.,6.},{6.,6.,6.}}),
    dealii::Tensor<2,3>({{6.,6.,6.},{6.,6.,6.},{6.,6.,6.}}),
  };
  for (unsigned int i = 0; i < result.size(); i++)
    {
      REQUIRE(result[i][0][0] == Approx(result_ref[i][0][0]));
    }
}
