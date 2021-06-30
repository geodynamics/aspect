/*
  Copyright (C) 2020 - 2021 by the authors of the ASPECT code.

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
#include <aspect/mesh_refinement/isosurfaces.h>

namespace aspect
{
  namespace MeshRefinement
  {
    namespace internal
    {
      // Forward declaration
      unsigned int min_max_string_to_int(const std::string &string_value, const unsigned int minimum_refinement_level, const unsigned int  maximum_refinement_level);
    }
  }
}

TEST_CASE("Isosurfaces min_max_string_to_int function")
{

  using namespace aspect::MeshRefinement::internal;
  CHECK(min_max_string_to_int("min",1,5) == 1);
  CHECK(min_max_string_to_int("min+1",1,5) == 2);
  CHECK(min_max_string_to_int("min+2",1,5) == 3);
  CHECK(min_max_string_to_int("min+3",1,5) == 4);
  CHECK(min_max_string_to_int("min+4",1,5) == 5);
  CHECK(min_max_string_to_int("max",1,5) == 5);
  CHECK(min_max_string_to_int("max-1",1,5) == 4);
  CHECK(min_max_string_to_int("max-2",1,5) == 3);
  CHECK(min_max_string_to_int("max-3",1,5) == 2);
  CHECK(min_max_string_to_int("max-4",1,5) == 1);
  CHECK(min_max_string_to_int("2",1,5) == 2);

  CHECK(min_max_string_to_int("6",1,5) == 5);
  CHECK(min_max_string_to_int("0",1,5) == 1);
  CHECK(min_max_string_to_int("min+5",1,5) == 5);
  CHECK(min_max_string_to_int("max-5",1,5) == 1);

  CHECK_THROWS_WITH(min_max_string_to_int("min-1",1,5),Contains("A value of min-1 was provided, but you can't provide a smaller value"));
  CHECK_THROWS_WITH(min_max_string_to_int("max+1",1,5),Contains("A value of max+1 was provided, but you can't provide a larger value"));
}
