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

TEST_CASE("floating point check example")
{
  // this would fail because of floating point arithmetic, but using Approx(f) will work:
  //  REQUIRE(0.1*3.0 == 0.3);
  REQUIRE(0.1*3.0 == Approx(0.3));
}

TEST_CASE("floating point std::vector checks")
{
  // here is how you can compare two std::vector<double>:
  std::vector<double> computed = {0.1*3.0, 1.0};
  std::vector<double> expected = {0.3, 1.0};

  compare_vectors_approx(computed, expected);
}

TEST_CASE("Utilities::read_and_distribute_file_content")
{
  const std::string content =
    aspect::Utilities::read_and_distribute_file_content(ASPECT_SOURCE_DIR "/VERSION",
                                                        MPI_COMM_WORLD);
  REQUIRE(content.size()>0);

  REQUIRE_THROWS_WITH(
    aspect::Utilities::read_and_distribute_file_content("does-not-exist.cc",
                                                        MPI_COMM_WORLD),
    Contains("Could not open"));
}
