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


#include <aspect/gravity_model/radial_linear.h>

#include <aspect/geometry_model/interface.h>
#include <aspect/coordinate_systems.h>

#include <deal.II/base/tensor.h>

namespace aspect
{
  namespace GravityModel
  {
    template <int dim>
    Tensor<1,dim>
    RadialLinear<dim>::gravity_vector (const Point<dim> &p) const
    {
      if (p.norm() == 0.0)
        return Tensor<1,dim>();

      const double depth = this->get_geometry_model().depth(p);
      const double maximal_depth = this->get_geometry_model().maximal_depth();

      return  (-magnitude_at_surface * (1.0 - depth/maximal_depth)
               -magnitude_at_bottom  * (depth/maximal_depth))
              * p/p.norm();
    }



    template <int dim>
    void
    RadialLinear<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Radial linear");
        {
          prm.declare_entry ("Magnitude at surface", "9.8",
                             Patterns::Double (),
                             "Magnitude of the radial gravity vector "
                             "at the surface of the domain. "
                             "Units: \\si{\\meter\\per\\second\\squared}.");
          prm.declare_entry ("Magnitude at bottom", "10.7",
                             Patterns::Double (),
                             "Magnitude of the radial gravity vector "
                             "at the bottom of the domain. `Bottom' means the "
                             "maximum depth in the chosen geometry, and for "
                             "example represents the core-mantle boundary in "
                             "the case of the `spherical shell' geometry model, "
                             "and the center in the case of the `sphere' "
                             "geometry model. Units: \\si{\\meter\\per\\second\\squared}.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }



    template <int dim>
    void
    RadialLinear<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Radial linear");
        {
          magnitude_at_surface = prm.get_double ("Magnitude at surface");
          magnitude_at_bottom = prm.get_double ("Magnitude at bottom");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      AssertThrow (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical ||
                   this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::ellipsoidal,
                   ExcMessage ("Gravity model 'radial linear' should not be used with geometry models that "
                               "do not have either a spherical or ellipsoidal natural coordinate system."));

    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace GravityModel
  {
    ASPECT_REGISTER_GRAVITY_MODEL(RadialLinear,
                                  "radial linear",
                                  "A gravity model which is radial (pointing inward if the gravity "
                                  "is positive) and the magnitude changes linearly with depth. The "
                                  "magnitude of gravity at the surface and bottom is read from the "
                                  "input file in a section ``Gravity model/Radial linear''.")
  }
}
