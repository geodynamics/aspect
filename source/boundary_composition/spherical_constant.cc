/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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


namespace aspect
{
  namespace BoundaryComposition
  {
    template <int dim>
    void
    SphericalConstant<dim>::
    initialize ()
    {
      // verify that the geometry is supported by this plugin
      AssertThrow ( Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>>(this->get_geometry_model()) ||
                    Plugins::plugin_type_matches<const GeometryModel::Sphere<dim>>(this->get_geometry_model()) ||
                    Plugins::plugin_type_matches<const GeometryModel::Chunk<dim>>(this->get_geometry_model()) ||
                    Plugins::plugin_type_matches<const GeometryModel::EllipsoidalChunk<dim>>(this->get_geometry_model()),
                    ExcMessage ("This boundary model is only implemented if the geometry is "
                                "one of the spherical geometries."));

      // no inner boundary in a full sphere
      if (Plugins::plugin_type_matches<const GeometryModel::Sphere<dim>>(this->get_geometry_model()))
        inner_boundary_indicator = numbers::invalid_unsigned_int;
      else
        inner_boundary_indicator = this->get_geometry_model().translate_symbolic_boundary_name_to_id("bottom");

      outer_boundary_indicator = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");

    }



    template <int dim>
    double
    SphericalConstant<dim>::
    boundary_composition (const types::boundary_id boundary_indicator,
                          const Point<dim> &/*position*/,
                          const unsigned int compositional_field) const
    {
      if (boundary_indicator == outer_boundary_indicator)
        return outer_composition[compositional_field];
      else if (boundary_indicator == inner_boundary_indicator)
        return inner_composition[compositional_field];
      else
        AssertThrow (false,
                     ExcMessage ("Unknown boundary indicator for geometry model. "
                                 "The given boundary should be ``top'' or ``bottom''."));

      return numbers::signaling_nan<double>();
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
                             Patterns::List(Patterns::Double ()),
                             "A comma separated list of composition boundary values "
                             "at the top boundary (at maximal radius). This list must have "
                             "one entry or as many entries as there are compositional fields. "
                             "Units: none.");
          prm.declare_entry ("Inner composition", "1.",
                             Patterns::List(Patterns::Double ()),
                             "A comma separated list of composition boundary values "
                             "at the bottom boundary (at minimal radius). This list must have "
                             "one entry or as many entries as there are compositional fields. "
                             "Units: none.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    SphericalConstant<dim>::parse_parameters (ParameterHandler &prm)
    {
      const unsigned int n_compositional_fields = this->n_compositional_fields();

      prm.enter_subsection("Boundary composition model");
      {
        prm.enter_subsection("Spherical constant");
        {
          inner_composition = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Inner composition"))),
                                                                      n_compositional_fields,
                                                                      "Inner boundary composition values");

          outer_composition = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Outer composition"))),
                                                                      n_compositional_fields,
                                                                      "Outer boundary composition values");
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
                                               "the inner and outer boundaries of a sphere, spherical "
                                               "shell, chunk or ellipsoidal chunk. "
                                               "Parameters are read from subsection 'Spherical constant'.")
  }
}
