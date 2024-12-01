/*
  Copyright (C) 2021 by the authors of the World Builder code.

  This file is part of the World Builder.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#define DOCTEST_CONFIG_SUPER_FAST_ASSERTS

#include "doctest/doctest.h"

#include "world_builder/utilities.h"

using namespace WorldBuilder;
using doctest::Approx;


TEST_CASE("multiply 3x3 matrices")
{
  std::array<std::array<double,3>,3> mat1 = {{{{1,2,3}},{{4,5,6}},{{7,8,9}}}};
  std::array<std::array<double,3>,3> mat2 = {{{{9,8,7}},{{6,5,4}},{{3,2,1}}}};
  std::array<std::array<double,3>,3> mat3 = {{{{0.1,10,-7}},{{4,-0.5,9}},{{-9,0.3,-5}}}};
  std::array<std::array<double,3>,3> result1 = {{{{30,24,18}},{{84,69,54}},{{138,114,90}}}}; //mat1*mat2
  std::array<std::array<double,3>,3> result2 = {{{{90,114,138}},{{54,69,84}},{{18,24,30}}}}; //mat2*mat1
  std::array<std::array<double,3>,3> result3 = {{{{-18.9,9.9,-4}},{{-33.6,39.3,-13}},{{-48.3,68.7,-22}}}}; //mat1*mat3
  std::array<std::array<double,3>,3> result4 = {{{{-8.9,-5.8,-2.7}},{{65,77.5,90}},{{-42.8,-56.5,-70.2}}}}; //mat3*mat1
  std::array<std::array<double,3>,3> result5 = {{{{-30.1,88.1,-26}},{{-15.4,58.7,-17}},{{-0.7,29.3,-8}}}}; //mat2*mat3
  std::array<std::array<double,3>,3> result6 = {{{{39.9,36.8,33.7}},{{60,47.5,35}},{{-94.2,-80.5,-66.8}}}}; //mat3*mat2

  std::array<std::array<double,3>,3> result7 = Utilities::multiply_3x3_matrices(mat1, mat2);
  for (size_t i = 0; i < 3; i++)
    {
      for (size_t j = 0; j < 3; j++)
        {
          CHECK(result1[i][j] == Approx(result7[i][j]));
        }
    }

  std::array<std::array<double,3>,3> result8 = Utilities::multiply_3x3_matrices(mat2, mat1);
  for (size_t i = 0; i < 3; i++)
    {
      for (size_t j = 0; j < 3; j++)
        {
          CHECK(result2[i][j] == Approx(result8[i][j]));
        }
    }

  std::array<std::array<double,3>,3> result9 = Utilities::multiply_3x3_matrices(mat1, mat3);
  for (size_t i = 0; i < 3; i++)
    {
      for (size_t j = 0; j < 3; j++)
        {
          CHECK(result3[i][j] == Approx(result9[i][j]));
        }
    }

  std::array<std::array<double,3>,3> result10 = Utilities::multiply_3x3_matrices(mat3, mat1);
  for (size_t i = 0; i < 3; i++)
    {
      for (size_t j = 0; j < 3; j++)
        {
          CHECK(result4[i][j] == Approx(result10[i][j]));
        }
    }

  std::array<std::array<double,3>,3> result11 = Utilities::multiply_3x3_matrices(mat2, mat3);
  for (size_t i = 0; i < 3; i++)
    {
      for (size_t j = 0; j < 3; j++)
        {
          CHECK(result5[i][j] == Approx(result11[i][j]));
        }
    }

  std::array<std::array<double,3>,3> result12 = Utilities::multiply_3x3_matrices(mat3, mat2);
  for (size_t i = 0; i < 3; i++)
    {
      for (size_t j = 0; j < 3; j++)
        {
          CHECK(result6[i][j] == Approx(result12[i][j]));
        }
    }

}