#include <aspect/gravity_model/vertical.h>

#include <deal.II/base/tensor.h>

namespace aspect
{
  namespace GravityModel
  {
    template <int dim>
    Tensor<1,dim>
    Vertical<dim>::gravity_vector (const Point<dim> &) const
    {
      Tensor<1,dim> g;
      g[dim-1] = -1;
      return g;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace GravityModel
  {
    template class Vertical<deal_II_dimension>;
    ASPECT_REGISTER_GRAVITY_MODEL(Vertical,
                                  "vertical",
                                  "A gravity model in which the gravity direction is vertically downward "
                                  "and at constant magnitude.");
  }
}
