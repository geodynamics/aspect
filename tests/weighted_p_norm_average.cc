//#include <aspect/simulator.h>
//#include <aspect/geometry_model/initial_topography_model/interface.h>
//#include <aspect/geometry_model/initial_topography_model/prm_polygon.h>
//#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

//#include <iostream>
//#include <typeinfo>

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

