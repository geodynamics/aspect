#include "common.h"
#include <aspect/utilities.h>

TEST_CASE("floating point checks")
{
  //std::cout << "hello!" <<std::endl;
  //  REQUIRE(0.1*3.0 == 0.3);
  REQUIRE(0.1*3.0 == Approx(0.3));
}

TEST_CASE("Utilities weighted_p_norm_average")
{
  std::vector<double> weights = {1,1,2,2,3,3};
  std::vector<double> values = {6,5,4,3,2,1};
  std::vector<double> p_norms = {-1000,-2.5,-2,-1,0,1,2,2.5,3,4,1000};
  std::vector<double> expected = {1., 1.59237, 1.6974 , 1.98895, 2.38899, 2.83333, 3.24037, 3.41824, 3.57872, 3.85347, 6. };
  
  for (unsigned int i = 0; i < p_norms.size(); i++)
    {
      INFO("check i=" << i << ": "); 
      //REQUIRE(aspect::Utilities::weighted_p_norm_average(weights,values,p_norms[i]) == expected[i]);
      REQUIRE(aspect::Utilities::weighted_p_norm_average(weights,values,p_norms[i]) == Approx(expected[i]));
    }
  
}
