// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2019-2024 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include "common.h"
#include <aspect/termination_criteria/steady_temperature.h>


/**
 * A function that trims the handed over list and removes all entries from the front that are
 * further back in time measured from the last entry than given by the first argument.
 * Additionally it makes sure to always keep two entries in the list, if the list had
 * two or more entries. Otherwise the function does not change the list.
 */
void trim_time_temperature_list (const double necessary_time_in_steady_state,
                                 std::list<std::pair<double, double>> &time_temperature_list)
{
  // Remove old times until we're at the correct time period
  // but ensure at least two entries remain in the list (one old, one current timestep)
  auto it = time_temperature_list.begin();
  while (time_temperature_list.back().first - (*it).first > necessary_time_in_steady_state &&
         std::distance(it,time_temperature_list.end()) > 2)
    ++it;

  time_temperature_list.erase(time_temperature_list.begin(), it);
}


TEST_CASE("trim time_temperature_list test for steady state termination criteria")
{
  std::list<std::pair<double,double>> time_temperature;
  time_temperature.emplace_back(0.0,1.0);
  time_temperature.emplace_back(1.0,2.0);
  time_temperature.emplace_back(1.4,3.0);

  trim_time_temperature_list(1.5,time_temperature);

  REQUIRE(time_temperature.front().second == 1.0);
  REQUIRE(time_temperature.back().second == 3.0);

  // Remove first entry because it is more than 1.0 away from last entry
  trim_time_temperature_list(1.0,time_temperature);

  REQUIRE(time_temperature.front().second == 2.0);
  REQUIRE(time_temperature.back().second == 3.0);

  // Try to remove first entry, should not change list, since it should keep at least 2 entries
  trim_time_temperature_list(0.3,time_temperature);

  REQUIRE(time_temperature.front().second == 2.0);
  REQUIRE(time_temperature.back().second == 3.0);
}
