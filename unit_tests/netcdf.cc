// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2021-2024 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

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
