#include <aspect/utilities.h>

int f()
{
  using namespace aspect;

  std::vector<double> weights = {1,2};
  std::vector<double> values = {3,4};
  std::vector<double> derivatives = {5,6};
  std::vector<double> p_norms = {-1000,-2.5,-2,-1,0,1,2,2.5,3,4,1000};

  for (unsigned int i = 0; i < p_norms.size(); i++)
    {
      std::cout << "p = " << p_norms[i] << ", average = " << aspect::Utilities::derivative_of_weighted_p_norm_average<double>(0,weights,values,derivatives,p_norms[i]) << std::endl;
    }

  exit(0);
  return 42;
}
// run this function by initializing a global variable by it
int i = f();

