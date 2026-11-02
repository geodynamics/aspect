// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2018-2022 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_unit_tests_common_h
#define _aspect_unit_tests_common_h

#include <catch.hpp>

#include <deal.II/base/table.h>


using Catch::Matchers::Contains;
using Catch::Matchers::StartsWith;

/**
 * Compare the given two std::vector<double> entries with an epsilon (using Catch::Approx)
 */
inline void compare_vectors_approx(
  const std::vector<double> &computed,
  const std::vector<double> &expected)
{
  REQUIRE(computed.size() == expected.size());
  for (unsigned int i=0; i<computed.size(); ++i)
    {
      INFO("array index i=" << i << ": ");
      REQUIRE(computed[i] == Approx(expected[i]));
    }
}

/**
 * Compare the given two table entries with an epsilon (using Catch::Approx)
 */
inline void compare_tables_approx(
  const dealii::Table<2,double> &computed,
  const dealii::Table<2,double> &expected)
{
  REQUIRE(computed.size() == expected.size());

  for (unsigned int i=0; i<computed.n_rows(); ++i)
    {
      for (unsigned int j=0; j<computed.n_cols(); ++j)
        {
          INFO("array index i,j=" << i << ',' << j << ": ");
          REQUIRE(computed[i][j] == Approx(expected[i][j]));
        }
    }
}

#endif
