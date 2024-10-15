/*
  Copyright (C) 2021 - 2024 by the authors of the ASPECT code.

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

TEST_CASE("Utilities::load_netcdf-1d")
{
  using namespace dealii;

  aspect::Utilities::StructuredDataLookup<1> lookup(1.0 /*scaling*/);
  lookup.load_netcdf(ASPECT_SOURCE_DIR "/data/test/netcdf/PREM_1s_depth.nc");

  REQUIRE(lookup.get_column_names().size()==9);
  REQUIRE(lookup.get_column_names()[0] == "radius");
  REQUIRE(lookup.get_column_names()[2] == "vpv");

  REQUIRE(lookup.get_data(Point<1>(1.0), 0) == Approx(6370.0)); // radius at depth = 1.0
  REQUIRE(lookup.get_data(Point<1>(1.0), 2) == Approx(1.45)); // vpv at depth = 1.0
}

TEST_CASE("Utilities::load_netcdf-1d-column-names")
{
  using namespace dealii;

  aspect::Utilities::StructuredDataLookup<1> lookup(1.0 /*scaling*/);
  lookup.load_netcdf(ASPECT_SOURCE_DIR "/data/test/netcdf/PREM_1s_depth.nc",
                     std::vector<std::string> {"vpv","radius"});

  REQUIRE(lookup.get_column_names().size()==2);
  REQUIRE(lookup.get_column_names()[0] == "vpv");
  REQUIRE(lookup.get_column_names()[1] == "radius");

  REQUIRE(lookup.get_data(Point<1>(1.0), 0) == Approx(1.45)); // vpv at depth = 1.0
  REQUIRE(lookup.get_data(Point<1>(1.0), 1) == Approx(6370.0)); // radius at depth = 1.0
}

TEST_CASE("Utilities::load_netcdf-3d")
{
  using namespace dealii;

  aspect::Utilities::StructuredDataLookup<3> lookup(1.0 /*scaling*/);
  lookup.load_netcdf(ASPECT_SOURCE_DIR "/data/test/netcdf/test-3d-cartesian.nc");

  REQUIRE(lookup.get_column_names().size()==2);
  REQUIRE(lookup.get_column_names()[0] == "index");
  REQUIRE(lookup.get_column_names()[1] == "depth");

  REQUIRE(lookup.get_data(Point<3>(0., 500., 0.), 0) == Approx(0.));
  REQUIRE(lookup.get_data(Point<3>(1000., 500., 0.), 0) == Approx(1.));
}

TEST_CASE("Utilities::load_netcdf-3d-spherical")
{
  using namespace dealii;

  aspect::Utilities::StructuredDataLookup<3> lookup(1.0 /*scaling*/);
  lookup.load_netcdf(ASPECT_SOURCE_DIR "/data/test/netcdf/test-3d-spherical.nc");

  REQUIRE(lookup.get_column_names().size()==1);
  REQUIRE(lookup.get_column_names()[0] == "vs_anomaly");

  REQUIRE(lookup.get_data(Point<3>(5297311.3326, 2.0944, 2.793), 0) == Approx(6423.5));
  REQUIRE(lookup.get_data(Point<3>(6252530.6549, 4.1888, 3.142), 0) == Approx(4574.6));
}

TEST_CASE("Utilities::load_netcdf-3d-column-names")
{
  using namespace dealii;

  aspect::Utilities::StructuredDataLookup<3> lookup(1.0 /*scaling*/);
  lookup.load_netcdf(ASPECT_SOURCE_DIR "/data/test/netcdf/test-3d-cartesian.nc",
                     std::vector<std::string> {"index"});

  REQUIRE(lookup.get_column_names().size()==1);
  REQUIRE(lookup.get_column_names()[0] == "index");

  REQUIRE(lookup.get_data(Point<3>(0., 500., 0.), 0) == Approx(0.));
  REQUIRE(lookup.get_data(Point<3>(1000., 500., 0.), 0) == Approx(1.));
}

#endif
