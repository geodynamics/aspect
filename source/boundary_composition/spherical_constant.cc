/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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
                          const unsigned int compositional_field) const
    {
      Assert(compositional_field == 0,
             ExcMessage("The 'spherical constant' boundary composition plugin "
                        "only supports a single compositional field."));
      (void) compositional_field;

      Assert (this->get_geometry_model().translate_id_to_symbol_name(boundary_indicator) == "top" ||
              this->get_geometry_model().translate_id_to_symbol_name(boundary_indicator) == "bottom",
              ExcMessage ("Unknown boundary indicator for geometry model. "
                          "The given boundary should be ``top'' or ``bottom''."));

      return boundary_compositions[boundary_indicator];
    }



    template <int dim>
    void
    SphericalConstant<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary composition model");
      {
        prm.enter_subsection("Spherical constant");
        {
          prm.declare_entry ("Outer composition", "0.",
                             Patterns::Double (),
                             "Composition at the outer boundary (lithosphere water/air). "
                             "For a spherical geometry model, this is the only boundary. "
                             "Units: none.");
          prm.declare_entry ("Inner composition", "1.",
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
          const auto boundary_names_map = this->get_geometry_model().get_symbolic_boundary_names_map();
          boundary_compositions.resize(boundary_names_map.size(),numbers::signaling_nan<double>());

          if (boundary_names_map.find("bottom") != boundary_names_map.end())
            boundary_compositions[boundary_names_map.at("bottom")] =  prm.get_double ("Inner composition");

          if (boundary_names_map.find("top") != boundary_names_map.end())
            boundary_compositions[boundary_names_map.at("top")] =  prm.get_double ("Outer composition");
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
