//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/gravity_model/radial.h>

#include <deal.II/base/tensor.h>

namespace aspect
{
  namespace GravityModel
  {
// ------------------------------ RadialConstant -------------------
    template <int dim>
    Tensor<1,dim>
    RadialConstant<dim>::gravity_vector (const Point<dim> &p) const
    {
      const double r = p.norm();
      return -magnitude * p/r;
    }



    template <int dim>
    void
    RadialConstant<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Radial constant");
        {
          prm.declare_entry ("Magnitude", "30",
                             Patterns::Double (0),
                             "Magnitude of the gravity vector in $m/s^2$. The direction is "
                             "always radially outward from the center of the earth.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    RadialConstant<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Radial constant");
        {
          magnitude = prm.get_double ("Magnitude");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }

// ------------------------------ RadialEarthLike -------------------

    template <int dim>
    Tensor<1,dim>
    RadialEarthLike<dim>::gravity_vector (const Point<dim> &p) const
    {
      const double r = p.norm();
      return -(1.245e-6 * r + 7.714e13/r/r) * p / r;
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace GravityModel
  {
    template class RadialConstant<deal_II_dimension>;
    ASPECT_REGISTER_GRAVITY_MODEL(RadialConstant,
                                  "radial constant",
                                  "A gravity model in which the gravity direction is radially inward "
                                  "and at constant magnitude. The magnitude is read from the parameter "
                                  "file in subsection 'Radial constant'.");

    template class RadialEarthLike<deal_II_dimension>;
    ASPECT_REGISTER_GRAVITY_MODEL(RadialEarthLike,
                                  "radial earth-like",
                                  "A gravity model in which the gravity direction is radially inward "
                                  "and with a magnitude that matches that of the earth at "
                                  "the core-mantle boundary as well as at the surface and "
                                  "in between is physically correct under the assumption "
                                  "of a constant density.");
  }
}
