// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2022-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

//
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
