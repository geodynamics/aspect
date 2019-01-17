/*
  Copyright (C) 2019 by the authors of the ASPECT code.

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
#include <aspect/termination_criteria/steady_temperature.h>

TEST_CASE("trim time_temperature_list test for steady state termination criteria")
{
  std::list<std::pair<double,double> > time_temperature;
  time_temperature.push_back(std::make_pair(0.0,1.0));
  time_temperature.push_back(std::make_pair(1.0,2.0));
  time_temperature.push_back(std::make_pair(1.4,3.0));

  aspect::TerminationCriteria::internal::trim_time_temperature_list(1.5,time_temperature);

  REQUIRE(time_temperature.front().second == 1.0);
  REQUIRE(time_temperature.back().second == 3.0);

  // Remove first entry because it is more than 1.0 away from last entry
  aspect::TerminationCriteria::internal::trim_time_temperature_list(1.0,time_temperature);

  REQUIRE(time_temperature.front().second == 2.0);
  REQUIRE(time_temperature.back().second == 3.0);

  // Try to remove first entry, should not change list, since it should keep at least 2 entries
  aspect::TerminationCriteria::internal::trim_time_temperature_list(0.3,time_temperature);

  REQUIRE(time_temperature.front().second == 2.0);
  REQUIRE(time_temperature.back().second == 3.0);
}
