#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/gravity_model/interface.h>

using namespace dealii;


template <int dim>
class MyGravity :
  public aspect::GravityModel::Interface<dim>
{
  public:
    virtual Tensor<1,dim> gravity_vector (const Point<dim> &position) const
      {
	Tensor<1,dim> ret;
	ret[0] = position[1];
	ret[1] = 42.0;
	return ret;
      }
};


// explicit instantiation
ASPECT_REGISTER_GRAVITY_MODEL(MyGravity, "my gravity", "no description")
