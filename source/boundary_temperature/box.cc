#include <aspect/boundary_temperature/box.h>
#include <aspect/geometry_model/box.h>

#include <utility>
#include <limits>


namespace aspect
{
  namespace BoundaryTemperature
  {
// ------------------------------ Box -------------------

    template <int dim>
    double
    Box<dim>::
    temperature (const GeometryModel::Interface<dim> &geometry_model,
                 const unsigned int                   boundary_indicator,
                 const Point<dim>                    &location) const
    {
      // verify that the geometry is in fact a box since only
      // for this geometry do we know for sure what boundary indicators it
      // uses and what they mean
      Assert (dynamic_cast<const GeometryModel::Box<dim>*>(&geometry_model)
              != 0,
              ExcMessage ("This boundary model is only implemented if the geometry is "
                          "in fact a box."));

      switch (boundary_indicator)
        {
          case 0:
            return 1;
          case 1:
            return 0;
          default:
            Assert (false, ExcMessage ("Unknown boundary indicator."));
            return std::numeric_limits<double>::quiet_NaN();
        }
    }


    template <int dim>
    double
    Box<dim>::
    minimal_temperature () const
    {
      return 0;
    }



    template <int dim>
    double
    Box<dim>::
    maximal_temperature () const
    {
      return 1;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryTemperature
  {
    template class Box<deal_II_dimension>;
    ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL(Box,
                                               "box",
                                               "A model in which the temperature is chosen constant on "
                                               "the left and right sides of a box.");
  }
}
