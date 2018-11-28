/*
  Copyright (C) 2018 by the authors of the ASPECT code.

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

TEST_CASE("Utilities::weighted_p_norm_average")
{
  std::vector<double> weights = {1,1,2,2,3,3};
  std::vector<double> values = {6,5,4,3,2,1};
  std::vector<double> p_norms = {-1000,-2.5,-2,-1,0,1,2,2.5,3,4,1000};
  std::vector<double> expected = {1., 1.59237, 1.6974 , 1.98895, 2.38899, 2.83333, 3.24037, 3.41824, 3.57872, 3.85347, 6. };

  for (unsigned int i = 0; i < p_norms.size(); i++)
    {
      INFO("check i=" << i << ": ");
      REQUIRE(aspect::Utilities::weighted_p_norm_average(weights,values,p_norms[i]) == Approx(expected[i]));
    }

}

//	1.) Test that the values coming back match what's expected
//	2.) Test that the stringstream is formated correctly ("POINTS", linebreaks, values ordered properly)
#if HAVE_LIBDAP
TEST_CASE("Utilities::read_and_distribute_file_content")
{
	MPI_Comm                            mpi_communicator;

	std::string expectedString = "# POINTS: 3 3 1 2 3 4 5 6 7 8 9";
	std::string returnString = aspect::Utilities::read_and_distribute_file_content("http://test.opendap.org/opendap/BALTO/test.balto.csv", mpi_communicator);
	returnString = returnString.replace(13, 1, "");
	returnString = returnString.replace(19, 1, " ");
	returnString = returnString.replace(25, 1, " ");
	returnString = returnString.replace(31, 1, "");
	INFO("check if the stringstream from read_and_distribute_file_content match the file text");
	REQUIRE( returnString == expectedString);
}
#endif
