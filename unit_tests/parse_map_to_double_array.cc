/*
  Copyright (C) 2018 - 2023 by the authors of the ASPECT code.

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


TEST_CASE("Utilities::MapParsing::parse_map_to_double_array")
{
  // Parse multicomponent properties
  {
    INFO("check 1: ");
    aspect::Utilities::MapParsing::Options options({"background","C1","C2","C3","C4","C5"}, "TestField");
    compare_vectors_approx(aspect::Utilities::MapParsing::parse_map_to_double_array ("C1:100, C2:200, C3:300, C4:400, C5:500, background:0",
                           options),
    {0.0,100.0,200.0,300.0,400.0,500.0});

    INFO("check 2: ");
    compare_vectors_approx(aspect::Utilities::MapParsing::parse_map_to_double_array ("all:100",options),
    {100.0,100.0,100.0,100.0,100.0,100.0});

    INFO("check 3: ");
    compare_vectors_approx(aspect::Utilities::MapParsing::parse_map_to_double_array ("background:0, C1:100, C2:200, C3:300, C4:400, C5:500",
                           options),
    {0.0,100.0,200.0,300.0,400.0,500.0});

    INFO("check 4: ");
    compare_vectors_approx(aspect::Utilities::MapParsing::parse_map_to_double_array ("background:0, C2:200, C1:100, C4:400, C5:500, C3:300",
                           options),
    {0.0,100.0,200.0,300.0,400.0,500.0});

    INFO("check 5: ");
    compare_vectors_approx(aspect::Utilities::MapParsing::parse_map_to_double_array ("C1:100, C2:200, background:0, C3:300, C4:400, C5:500",
                           options),
    {0.0,100.0,200.0,300.0,400.0,500.0});
  }

  {
    INFO("check 6: ");
    aspect::Utilities::MapParsing::Options options({"C1","C2","C3","C4","C5"}, "TestField");
    compare_vectors_approx(aspect::Utilities::MapParsing::parse_map_to_double_array ("100, 200, 300, 400, 500",
                           options),
    {100.0,200.0,300.0,400.0,500.0});

    INFO("check 7: ");
    compare_vectors_approx(aspect::Utilities::MapParsing::parse_map_to_double_array ("C1:100, C2:200, C3:300, C4:400, C5:500",
                           options),
    {100.0,200.0,300.0,400.0,500.0});
  }

  {
    INFO("check 8: ");
    aspect::Utilities::MapParsing::Options options({"background","C1","C2","C3","C4","C5"}, "TestField");
    compare_vectors_approx(aspect::Utilities::MapParsing::parse_map_to_double_array ("background:0, C2:200, C1:100, C4:400, C5:500, C3:300",
                           options),
    {0.0,100.0,200.0,300.0,400.0,500.0});
  }

  {
    INFO("check 9: ");
    aspect::Utilities::MapParsing::Options options({"background","C1","C2","C3","C4","C5"}, "TestField");
    options.store_values_per_key = true;
    options.allow_multiple_values_per_key = true;
    compare_vectors_approx(aspect::Utilities::MapParsing::parse_map_to_double_array ("C1:100, C2:200|100, C3:300, C4:400, C5:500, background:0",
                           options),
    {0.0,100.0,200.0,100.0,300.0,400.0,500.0});

    REQUIRE(options.n_values_per_key == std::vector<unsigned int>({1,1,2,1,1,1}));
  }

  {
    INFO("check 10: ");
    aspect::Utilities::MapParsing::Options options({"C1","C2","C3","C4","C5"}, "TestField");
    options.allow_multiple_values_per_key = true;
    options.store_values_per_key = true;
    compare_vectors_approx(aspect::Utilities::MapParsing::parse_map_to_double_array ("C1:100, C2:200|, C3:300, C4:400, C5:500",
                           options),
    {100.0,200.0,300.0,400.0,500.0});

    REQUIRE(options.n_values_per_key == std::vector<unsigned int>({1,1,1,1,1}));
  }

  {
    INFO("check 11: ");
    aspect::Utilities::MapParsing::Options options({"C1","C2","C3","C4","C5"}, "TestField");
    options.store_values_per_key = true;
    options.allow_multiple_values_per_key = true;

    compare_vectors_approx(aspect::Utilities::MapParsing::parse_map_to_double_array ("C1:100, C2:200|300, C3:300, C4:400, C5:500",
                           options),
    {100.0,200.0,300.0,300.0,400.0,500.0});

    REQUIRE(options.n_values_per_key == std::vector<unsigned int>({1,2,1,1,1}));

    options.store_values_per_key = false;
    options.check_values_per_key = true;
    compare_vectors_approx(aspect::Utilities::MapParsing::parse_map_to_double_array ("C1:100, C2:200|300, C3:300, C4:400, C5:500",
                           options),
    {100.0,200.0,300.0,300.0,400.0,500.0});
  }

  {
    INFO("check 12: ");
    aspect::Utilities::MapParsing::Options options({"C1","C2"}, "TestField");
    options.allow_multiple_values_per_key = true;
    options.store_values_per_key = true;

    compare_vectors_approx(aspect::Utilities::MapParsing::parse_map_to_double_array ("all:300|400",
                           options),
    {300.0,400.0,300.0,400.0});

    REQUIRE(options.n_values_per_key == std::vector<unsigned int> ({2,2}));

    options.store_values_per_key = false;
    options.check_values_per_key = true;
    compare_vectors_approx(aspect::Utilities::MapParsing::parse_map_to_double_array ("all:100|200",
                           options),
    {100.0,200.0,100.0,200.0});
  }

  {
    INFO("check 13: ");
    aspect::Utilities::MapParsing::Options options({"background","C1","C2","C3","C4","C5"}, "TestField");
    options.allow_missing_keys = true;

    compare_vectors_approx(aspect::Utilities::MapParsing::parse_map_to_double_array ("background:0, C1:100, C4:400, C5:500, C3:300",
                           options),
    {0.0,100.0,300.0,400.0,500.0});
  }

  {
    INFO("check 14: ");
    aspect::Utilities::MapParsing::Options options({"background","C1","C2","C3","C4","C5"}, "TestField");
    options.allow_missing_keys = true;

    compare_vectors_approx(aspect::Utilities::MapParsing::parse_map_to_double_array ("C1:100, C4:400, C5:500, C3:300",
                           options),
    {100.0,300.0,400.0,500.0});
  }

  {
    INFO("check 16: ");
    aspect::Utilities::MapParsing::Options options({"background","C1","C2","C3","C4","C5"}, "TestField");
    options.allow_missing_keys = true;

    compare_vectors_approx(aspect::Utilities::MapParsing::parse_map_to_double_array ("",
                           options),
                           std::vector<double>());
  }

  {
    INFO("check 17: ");
    aspect::Utilities::MapParsing::Options options({"C1","C2","C3","C4","C5"}, "TestField");
    options.allow_multiple_values_per_key = true;
    options.store_values_per_key = true;
    options.allow_missing_keys = true;

    compare_vectors_approx(aspect::Utilities::MapParsing::parse_map_to_double_array ("C1:100, C2:200|300, C3:300, C5:500",
                           options),
    {100.0,200.0,300.0,300.0,500.0});

    REQUIRE(options.n_values_per_key == std::vector<unsigned int>({1,2,1,0,1}));

    options.store_values_per_key = false;
    options.check_values_per_key = true;
    compare_vectors_approx(aspect::Utilities::MapParsing::parse_map_to_double_array ("C1:100, C2:200|300, C3:300, C5:500",
                           options),
    {100.0,200.0,300.0,300.0,500.0});
  }

  {
    INFO("check 18: ");

    aspect::Utilities::MapParsing::Options options({"C1","C2"}, "TestField");
    options.allow_multiple_values_per_key = true;
    options.check_values_per_key = true;
    options.n_values_per_key = {2,2};

    compare_vectors_approx(aspect::Utilities::MapParsing::parse_map_to_double_array ("C1:300|400, C2:200",
                           options),
    {300.0,400.0,200.0,200.0});

    INFO("check 19: ");

    // still using the compare_vectors_approx defined before
    compare_vectors_approx(aspect::Utilities::MapParsing::parse_map_to_double_array ("200",
                           options),
    {200.0,200.0,200.0,200.0});
  }

  {
    // Check that we can provide additional allowed keys that are ignored when parsing
    INFO("check 20: ");
    aspect::Utilities::MapParsing::Options options({"C1","C2","C3","C5"},"TestField");
    options.list_of_allowed_keys = {"C1","C2","C3","C4","C5","C6","C7"};

    std::string input_string = "C1:100, C2:200, C3:300, C5:500";

    compare_vectors_approx(aspect::Utilities::MapParsing::parse_map_to_double_array (input_string, options),
    {100.0,200.0,300.0,500.0});
  }
  INFO("check complete");

}

using Catch::Matchers::Contains;

TEST_CASE("Utilities::MapParsing::parse_map_to_double_array FAIL ON PURPOSE")
{
  // Purposefully fail Parse multicomponent properties
  {
    INFO("check fail 1: ");
    aspect::Utilities::MapParsing::Options options({"C1","C2","C3","C4","C5"}, "TestField");

    REQUIRE_THROWS_WITH(
      aspect::Utilities::MapParsing::parse_map_to_double_array ("C1:100, C1:200, C3:300, C4:400, C5:500",
                                                                options), Contains("is listed multiple times"));
  }

  {
    INFO("check fail 2: ");
    aspect::Utilities::MapParsing::Options options({"background","C1","C2","C3","C4","C5"}, "TestField");

    REQUIRE_THROWS_WITH(
      aspect::Utilities::MapParsing::parse_map_to_double_array ("C1:100, C2:200, C3:300, C4:400, C5:500",
                                                                options), Contains("The keyword <background> in TestField is not listed"));
  }

  {
    INFO("check fail 3: ");
    aspect::Utilities::MapParsing::Options options({"bg","C1","C2","C3","C4","C5"}, "TestField");

    REQUIRE_THROWS_WITH(
      aspect::Utilities::MapParsing::parse_map_to_double_array ("C1:100, C1:200, C3:300, C4;400, C5:500, bg:3",
                                                                options), Contains("does not have the expected format"));
  }

  {
    INFO("check fail 4: ");
    aspect::Utilities::MapParsing::Options options({"bg","C1","C2","C3","C4","C5"}, "TestField");

    REQUIRE_THROWS_WITH(
      aspect::Utilities::MapParsing::parse_map_to_double_array ("C1:100, C2:200, C3:300, C4:400, C24:500, bg:3",
                                                                options), Contains("does not match any entries"));
  }

  {
    INFO("check fail 5: ");
    aspect::Utilities::MapParsing::Options options({"bg","C1","C2","C3","C4","C5"}, "TestField");

    REQUIRE_THROWS_WITH(
      aspect::Utilities::MapParsing::parse_map_to_double_array ("C1:100, C2:200, C3:300, all:400, C5:500, bg:3",
                                                                options), Contains("The keyword `all' in the property TestField is only allowed if"));
  }

  {
    INFO("check fail 6: ");
    aspect::Utilities::MapParsing::Options options({"background","C1","C2","C3","C4","C5"}, "TestField");

    REQUIRE_THROWS_WITH(
      aspect::Utilities::MapParsing::parse_map_to_double_array ("C1:100",
                                                                options), Contains("The keyword <background> in TestField is not listed"));
  }

  {
    // Subentries not allowed
    INFO("check fail 7: ");
    aspect::Utilities::MapParsing::Options options({"C1","C2","C3","C4","C5"}, "TestField");

    REQUIRE_THROWS_WITH(
      aspect::Utilities::MapParsing::parse_map_to_double_array ("C1:100, C2:100|200, C3:300, C4:400, C5:500",
                                                                options), Contains("The keyword <C2> in TestField has multiple values"));
  }

  {
    // Wrong subentry format
    INFO("check fail 8: ");
    aspect::Utilities::MapParsing::Options options({"C1","C2","C3","C4","C5"}, "TestField");
    options.allow_multiple_values_per_key = true;

    REQUIRE_THROWS_WITH(
      aspect::Utilities::MapParsing::parse_map_to_double_array ("C1:100, C2:|200, C3:300, C4:400, C5:500",
                                                                options), Contains("does not have the expected format"));
  }

  {
    // No subentries
    INFO("check fail 9: ");
    aspect::Utilities::MapParsing::Options options({"C1","C2","C3","C4","C5"}, "TestField");
    options.allow_multiple_values_per_key = true;

    REQUIRE_THROWS_WITH(
      aspect::Utilities::MapParsing::parse_map_to_double_array ("C1:100, C2:|, C3:300, C4:400, C5:500",
                                                                options), Contains("does not have the expected format"));
  }

  {
    // Wrong input structure
    INFO("check fail 10: ");
    aspect::Utilities::MapParsing::Options options({"C1","C2","C3","C4","C5"}, "TestField");
    options.allow_multiple_values_per_key = true;
    options.store_values_per_key = true;

    compare_vectors_approx(aspect::Utilities::MapParsing::parse_map_to_double_array ("C1:100, C2:200|300, C3:300, C4:400, C5:500",
                           options),
    {100.0,200.0,300.0,300.0,400.0,500.0});

    REQUIRE(options.n_values_per_key == std::vector<unsigned int>({1,2,1,1,1}));

    options.store_values_per_key = false;
    options.check_values_per_key = true;

    REQUIRE_THROWS_WITH(aspect::Utilities::MapParsing::parse_map_to_double_array ("C1:100|200, C2:300, C3:300, C4:400, C5:500",
                                                                                  options),
                        Contains("The key <C1> in <TestField> does not have the expected number"));
  }

  {
    INFO("check fail 11: ");
    aspect::Utilities::MapParsing::Options options({"C1","C2"}, "TestField");
    options.allow_multiple_values_per_key = true;
    options.store_values_per_key = true;

    compare_vectors_approx(aspect::Utilities::MapParsing::parse_map_to_double_array ("all:300|400",
                           options),
    {300.0,400.0,300.0,400.0});

    REQUIRE(options.n_values_per_key == std::vector<unsigned int>({2,2}));

    options.store_values_per_key = false;
    options.check_values_per_key = true;

    REQUIRE_THROWS_WITH(aspect::Utilities::MapParsing::parse_map_to_double_array ("all:100|200|300",
                                                                                  options),
                        Contains("The key <C1> in <TestField> does not have the expected number"));
  }
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
