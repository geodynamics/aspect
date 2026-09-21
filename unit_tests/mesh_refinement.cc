// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2020-2021 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

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
