/*
  Copyright (C) 2018 by the authors of the ASPECT code.

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
#include <deal.II/base/exceptions.h>


TEST_CASE("Utilities::parse_map_to_double_array")
{
  // Parse multicomponent properties
  INFO("check 1: ");
  compare_vectors_approx(aspect::Utilities::parse_map_to_double_array ("C1:100, C2:200, C3:300, C4:400, C5:500, background:0",
  {"C1","C2","C3","C4","C5"},
  true,
  "TestField"),
  {0.0,100.0,200.0,300.0,400.0,500.0});

  INFO("check 2: ");
  compare_vectors_approx(aspect::Utilities::parse_map_to_double_array ("all:100",
  {"C1","C2","C3","C4","C5"},
  true,
  "TestField"),
  {100.0,100.0,100.0,100.0,100.0,100.0});

  INFO("check 3: ");
  compare_vectors_approx(aspect::Utilities::parse_map_to_double_array ("background:0, C1:100, C2:200, C3:300, C4:400, C5:500",
  {"C1","C2","C3","C4","C5"},
  true,
  "TestField"),
  {0.0,100.0,200.0,300.0,400.0,500.0});

  INFO("check 4: ");
  compare_vectors_approx(aspect::Utilities::parse_map_to_double_array ("background:0, C2:200, C1:100, C4:400, C5:500, C3:300",
  {"C1","C2","C3","C4","C5"},
  true,
  "TestField"),
  {0.0,100.0,200.0,300.0,400.0,500.0});

  INFO("check 5: ");
  compare_vectors_approx(aspect::Utilities::parse_map_to_double_array ("C1:100, C2:200, background:0, C3:300, C4:400, C5:500",
  {"C1","C2","C3","C4","C5"},
  true,
  "TestField"),
  {0.0,100.0,200.0,300.0,400.0,500.0});

  INFO("check 6: ");
  compare_vectors_approx(aspect::Utilities::parse_map_to_double_array ("100, 200, 300, 400, 500",
  {"C1","C2","C3","C4","C5"},
  false,
  "TestField"),
  {100.0,200.0,300.0,400.0,500.0});

  INFO("check 7: ");
  compare_vectors_approx(aspect::Utilities::parse_map_to_double_array ("C1:100, C2:200, C3:300, C4:400, C5:500",
  {"C1","C2","C3","C4","C5"},
  false,
  "TestField"),
  {100.0,200.0,300.0,400.0,500.0});

  INFO("check 8: ");
  compare_vectors_approx(aspect::Utilities::parse_map_to_double_array ("background:0, C2:200, C1:100, C4:400, C5:500, C3:300",
  {"C1","C2","C3","C4","C5"},
  true,
  "TestField"),
  {0.0,100.0,200.0,300.0,400.0,500.0});

  {
    INFO("check 9: ");
    auto subentry_structure = std::make_shared<std::vector<unsigned int>>();

    compare_vectors_approx(aspect::Utilities::parse_map_to_double_array ("C1:100, C2:200|100, C3:300, C4:400, C5:500, background:0",
    {"C1","C2","C3","C4","C5"},
    true,
    "TestField",
    true,
    subentry_structure),
    {0.0,100.0,200.0,100.0,300.0,400.0,500.0});

    REQUIRE(*subentry_structure == std::vector<unsigned int>({1,1,2,1,1,1}));
  }

  {
    INFO("check 10: ");
    auto subentry_structure = std::make_shared<std::vector<unsigned int>>();

    compare_vectors_approx(aspect::Utilities::parse_map_to_double_array ("C1:100, C2:200|, C3:300, C4:400, C5:500",
    {"C1","C2","C3","C4","C5"},
    false,
    "TestField",
    true,
    subentry_structure),
    {100.0,200.0,300.0,400.0,500.0});

    REQUIRE(*subentry_structure == std::vector<unsigned int>({1,1,1,1,1}));
  }

  INFO("check complete");

}

using Catch::Matchers::Contains;

TEST_CASE("Utilities::parse_map_to_double_array FAIL ON PURPOSE")
{
  // Purposefully fail Parse multicomponent properties
  INFO("check fail 1: ");
  REQUIRE_THROWS_WITH(
    aspect::Utilities::parse_map_to_double_array ("C1:100, C1:200, C3:300, C4:400, C5:500",
  {"C1","C2","C3","C4","C5"},
  false,
  "TestField"), Contains("is listed multiple times"));


  INFO("check fail 2: ");
  REQUIRE_THROWS_WITH(
    aspect::Utilities::parse_map_to_double_array ("C1:100, C2:200, C3:300, C4:400, C5:500",
  {"C1","C2","C3","C4","C5"},
  true,
  "TestField"), Contains("list must equal one of the following values:"));


  INFO("check fail 3: ");
  REQUIRE_THROWS_WITH(
    aspect::Utilities::parse_map_to_double_array ("C1:100, C1:200, C3:300, C4;400, C5:500, bg:3",
  {"C1","C2","C3","C4","C5"},
  true,
  "TestField"), Contains("The required format for field"));


  INFO("check fail 4: ");
  REQUIRE_THROWS_WITH(
    aspect::Utilities::parse_map_to_double_array ("C1:100, C2:200, C3:300, C4:400, C24:500, bg:3",
  {"C1","C2","C3","C4","C5"},
  true,
  "TestField"), Contains("does not match any entries"));

  INFO("check fail 5: ");
  REQUIRE_THROWS_WITH(
    aspect::Utilities::parse_map_to_double_array ("C1:100, C2:200, C3:300, all:400, C5:500, bg:3",
  {"C1","C2","C3","C4","C5"},
  true,
  "TestField"), Contains("values found, the keyword `all' is not"));

  INFO("check fail 6: ");
  REQUIRE_THROWS_WITH(
    aspect::Utilities::parse_map_to_double_array ("C1:100",
  {"C1","C2","C3","C4","C5"},
  true,
  "TestField"), Contains("The keyword `all' is expected but is not found"));

  // Subentries not allowed
  INFO("check fail 7: ");
  REQUIRE_THROWS_WITH(
    aspect::Utilities::parse_map_to_double_array ("C1:100, C2:|200, C3:300, C4:400, C5:500",
  {"C1","C2","C3","C4","C5"},
  false,
  "TestField"), Contains("The required format for field"));

  // Wrong subentry format
  INFO("check fail 7: ");
  REQUIRE_THROWS_WITH(
    aspect::Utilities::parse_map_to_double_array ("C1:100, C2:|200, C3:300, C4:400, C5:500",
  {"C1","C2","C3","C4","C5"},
  false,
  "TestField",
  true), Contains("The required format for field"));

  // No subentries
  INFO("check fail 8: ");
  REQUIRE_THROWS_WITH(
    aspect::Utilities::parse_map_to_double_array ("C1:100, C2:|, C3:300, C4:400, C5:500",
  {"C1","C2","C3","C4","C5"},
  false,
  "TestField",
  true), Contains("The required format for field"));
}


TEST_CASE("Utilities::parse_input_table")
{
  // Parse multicomponent properties
  INFO("check 1: ");

  dealii::Table<2,double> expected_result(3,3);

  expected_result[0][0] = 100.0;
  expected_result[0][1] = 10.0;
  expected_result[0][2] = 30.0;
  expected_result[1][0] = 25.0;
  expected_result[1][1] = 5.0;
  expected_result[1][2] = 1500.0;
  expected_result[2][0] = 0.001;
  expected_result[2][1] = 0.1;
  expected_result[2][2] = 3.5;


  compare_tables_approx(aspect::Utilities::parse_input_table<double> ("100,10,30; 25.0,5.0 ,1.5e3; 0.001,10e-2,3.5",
                                                                      3, 3, "TestField"),
                        expected_result);

  expected_result[0][0] = 100.0;
  expected_result[0][1] = 10.0;
  expected_result[0][2] = 30.0;
  expected_result[1][0] = 100.0;
  expected_result[1][1] = 10.0;
  expected_result[1][2] = 30.0;
  expected_result[2][0] = 100.0;
  expected_result[2][1] = 10.0;
  expected_result[2][2] = 30.0;


  compare_tables_approx(aspect::Utilities::parse_input_table<double> ("100,10,30",
                                                                      3, 3, "TestField"),
                        expected_result);

  expected_result[0][0] = 100.0;
  expected_result[0][1] = 100.0;
  expected_result[0][2] = 100.0;
  expected_result[1][0] = 10.0;
  expected_result[1][1] = 10.0;
  expected_result[1][2] = 10.0;
  expected_result[2][0] = 30.0;
  expected_result[2][1] = 30.0;
  expected_result[2][2] = 30.0;


  compare_tables_approx(aspect::Utilities::parse_input_table<double> ("100;10;30",
                                                                      3, 3, "TestField"),
                        expected_result);

  INFO("check complete");

}

using Catch::Matchers::Contains;

TEST_CASE("Utilities::parse_input_table FAIL ON PURPOSE")
{
  // Purposefully fail Parse multicomponent properties
  INFO("check fail 1: ");
  REQUIRE_THROWS_WITH(
    aspect::Utilities::parse_input_table<double> ("100,10,30; 25.0,5.0 ,1.5e3",
                                                  3,3,
                                                  "TestField"), Contains("Length of TestField"));

  // Purposefully fail Parse multicomponent properties
  INFO("check fail 2: ");
  REQUIRE_THROWS_WITH(
    aspect::Utilities::parse_input_table<double> ("100,10; 25.0,5.0; 5.0,2.5",
                                                  3,3,
                                                  "TestField"), Contains("Length of TestField"));

  // Purposefully fail Parse multicomponent properties
  INFO("check fail 3: ");
  REQUIRE_THROWS_WITH(
    aspect::Utilities::parse_input_table<double> ("100,10; a,5.0; 5.0,2.5",
                                                  3,2,
                                                  "TestField"), Contains("bad lexical cast"));
}
