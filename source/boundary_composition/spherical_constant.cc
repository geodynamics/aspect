/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/boundary_composition/spherical_constant.h>
#include <aspect/geometry_model/spherical_shell.h>

#include <utility>
#include <limits>


namespace aspect
{
  namespace BoundaryComposition
  {
// ------------------------------ SphericalConstant -------------------

    template <int dim>
    double
    SphericalConstant<dim>::
    composition (const GeometryModel::Interface<dim> &geometry_model,
                 const unsigned int                   boundary_indicator,
                 const Point<dim>                    &location,
                 const unsigned int                   compositional_field) const
    {
      // verify that the geometry is in fact a spherical shell since only
      // for this geometry do we know for sure what boundary indicators it
      // uses and what they mean
      Assert (dynamic_cast<const GeometryModel::SphericalShell<dim>*>(&geometry_model)
              != 0,
              ExcMessage ("This boundary model is only implemented if the geometry is "
                          "in fact a spherical shell."));

      switch (boundary_indicator)
        {
          case 0:
            return inner_composition;
          case 1:
            return outer_composition;
          default:
            Assert (false, ExcMessage ("Unknown boundary indicator."));
            return std::numeric_limits<double>::quiet_NaN();
        }
    }


    template <int dim>
    double
    SphericalConstant<dim>::
    minimal_composition (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      return std::min (inner_composition, outer_composition);
    }



    template <int dim>
    double
    SphericalConstant<dim>::
    maximal_composition (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      return std::max (inner_composition, outer_composition);
    }



    template <int dim>
    void
    SphericalConstant<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary composition model");
      {
        prm.enter_subsection("Spherical constant");
        {
          prm.declare_entry ("Outer composition", "0",
                             Patterns::Double (),
                             "Composition at the outer boundary (lithosphere water/air). Units: K.");
          prm.declare_entry ("Inner composition", "1",
                             Patterns::Double (),
                             "Composition at the inner boundary (core mantle boundary). Units: K.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    SphericalConstant<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary composition model");
      {
        prm.enter_subsection("Spherical constant");
        {
          inner_composition = prm.get_double ("Inner composition");
          outer_composition = prm.get_double ("Outer composition");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryComposition
  {
    ASPECT_REGISTER_BOUNDARY_COMPOSITION_MODEL(SphericalConstant,
                                               "spherical constant",
                                               "A model in which the composition is chosen constant on "
                                               "the inner and outer boundaries of a spherical shell. "
                                               "Parameters are read from subsection 'Sherical constant'.")
  }
}
