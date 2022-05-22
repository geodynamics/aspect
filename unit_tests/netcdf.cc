/*
  Copyright (C) 2021 - 2022 by the authors of the ASPECT code.

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

#ifdef ASPECT_WITH_NETCDF

TEST_CASE("Utilities::load_netcdf")
{
  using namespace dealii;

  aspect::Utilities::StructuredDataLookup<1> lookup(1.0 /*scaling*/);
  lookup.load_netcdf(ASPECT_SOURCE_DIR "/data/test/netcdf/PREM_1s_depth.nc");

  REQUIRE(lookup.get_column_names()[0] == "radius");
  REQUIRE(lookup.get_column_names()[2] == "vpv");

  REQUIRE(lookup.get_data(Point<1>(1.0), 0) == Approx(6370.0)); // radius at depth = 1.0
  REQUIRE(lookup.get_data(Point<1>(1.0), 2) == Approx(1.45)); // vpv at depth = 1.0
}


#endif
