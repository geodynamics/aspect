/*
  Copyright (C) 2022 by the authors of the ASPECT code.

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

//#include <aspect/simulator.h>
//#include <aspect/geometry_model/initial_topography_model/interface.h>
//#include <aspect/geometry_model/initial_topography_model/prm_polygon.h>
//#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

//#include <iostream>
//
int f()
{
  using namespace aspect;

  std::vector<double> weights = {1,1,2,2,3,3};
  std::vector<double> values = {6,5,4,3,2,1};
  std::vector<double> p_norms = {-1000,-2.5,-2,-1,0,1,2,2.5,3,4,1000};

  for (unsigned int i = 0; i < p_norms.size(); i++)
    {
      std::cout << "p = " << p_norms[i] << ", average = " << aspect::Utilities::weighted_p_norm_average(weights,values,p_norms[i]) << std::endl;
    }

  exit(0);
  return 42;
}
// run this function by initializing a global variable by it
int i = f();
