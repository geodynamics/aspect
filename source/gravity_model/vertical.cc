// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/gravity_model/vertical.h>

#include <aspect/geometry_model/interface.h>

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
      g[dim-1] = -gravity_magnitude;
      return g;
    }



    template <int dim>
    void
    Vertical<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Vertical");
        {
          prm.declare_entry ("Magnitude", "1.",
                             Patterns::Double (),
                             "Value of the gravity vector in $\\si{\\meter\\per\\second\\squared}$ directed "
                             "along the negative $y$ (in 2d) or $z$ (in 3d) axis (if the magnitude "
                             "is positive.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }



    template <int dim>
    void
    Vertical<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Vertical");
        {
          gravity_magnitude = prm.get_double ("Magnitude");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      AssertThrow (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::cartesian,
                   ExcMessage ("Gravity model 'vertical' should not be used with geometry models that "
                               "do not have a Cartesian natural coordinate system."));
    }
  }
}



// explicit instantiations
namespace aspect
{
  namespace GravityModel
  {
    ASPECT_REGISTER_GRAVITY_MODEL(Vertical,
                                  "vertical",
                                  "A gravity model in which the gravity direction is vertical (pointing "
                                  "downward for positive values) and at a constant magnitude by default "
                                  "equal to one.")
  }
}
