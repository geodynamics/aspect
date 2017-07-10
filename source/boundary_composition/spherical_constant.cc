/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <deal.II/base/signaling_nan.h>
#include <aspect/boundary_composition/spherical_constant.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/chunk.h>
#include <aspect/geometry_model/ellipsoidal_chunk.h>

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
    boundary_composition (const types::boundary_id boundary_indicator,
                          const Point<dim> &/*position*/,
                          const unsigned int /*compositional_field*/) const
    {
      const GeometryModel::Interface<dim> *geometry_model = &this->get_geometry_model();
      const std::string boundary_name = geometry_model->translate_id_to_symbol_name(boundary_indicator);

      if (boundary_name == "bottom")
        return inner_composition;
      else if (boundary_name =="top")
        return outer_composition;
      else
        {
          Assert (false, ExcMessage ("Unknown boundary indicator for geometry model. "
                                     "The given boundary should be ``top'' or ``bottom''."));
          return numbers::signaling_nan<double>();
        }
    }



    template <int dim>
    double
    SphericalConstant<dim>::
    minimal_composition (const std::set<types::boundary_id> &) const
    {
      const GeometryModel::Interface<dim> *geometry_model = &this->get_geometry_model();
      if (dynamic_cast<const GeometryModel::Sphere<dim>*>(geometry_model) != 0)
        return outer_composition;
      else
        return std::min (inner_composition, outer_composition);
    }



    template <int dim>
    double
    SphericalConstant<dim>::
    maximal_composition (const std::set<types::boundary_id> &) const
    {
      const GeometryModel::Interface<dim> *geometry_model = &this->get_geometry_model();
      if (dynamic_cast<const GeometryModel::Sphere<dim>*>(geometry_model) != 0)
        return outer_composition;
      else
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
                             "Composition at the outer boundary (lithosphere water/air). "
                             "For a spherical geometry model, this is the only boundary. "
                             "Units: none.");
          prm.declare_entry ("Inner composition", "1",
                             Patterns::Double (),
                             "Composition at the inner boundary (core mantle boundary). Units: none.");
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
                                               "the inner and outer boundaries of a surface, spherical "
                                               "shell, chunk or ellipsoidal chunk. "
                                               "Parameters are read from subsection 'Spherical constant'.")
  }
}
