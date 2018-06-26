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
  compare_vectors_approx(aspect::Utilities::parse_map_to_double_array ("C1:100, C2:200, C3:300, C4:400, C5:500, bg:0",
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
  compare_vectors_approx(aspect::Utilities::parse_map_to_double_array ("bg:0, C1:100, C2:200, C3:300, C4:400, C5:500",
  {"C1","C2","C3","C4","C5"},
  true,
  "TestField"),
  {0.0,100.0,200.0,300.0,400.0,500.0});

  INFO("check 4: ");
  compare_vectors_approx(aspect::Utilities::parse_map_to_double_array ("bg:0, C2:200, C1:100, C4:400, C5:500, C3:300",
  {"C1","C2","C3","C4","C5"},
  true,
  "TestField"),
  {0.0,100.0,200.0,300.0,400.0,500.0});

  INFO("check 5: ");
  compare_vectors_approx(aspect::Utilities::parse_map_to_double_array ("C1:100, C2:200, bg:0, C3:300, C4:400, C5:500",
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


  /*
    try
    {
        aspect::Utilities::parse_map_to_double_array ("C1:100, C1:200, C3:300, C4:400, C5:500",
    {"C1","C2","C3","C4","C5"},
    false,
    "TestField");
    }
    catch (dealii::ExceptionBase &e)
    {
      e.print_info(std::cerr);
      std::cerr << "waht=" << e.what() << std::endl;
      //INFO(e.what());
      REQUIRE(false);
    }
  */

}
