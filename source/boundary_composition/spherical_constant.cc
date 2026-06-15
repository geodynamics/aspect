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
    }



    template <int dim>
    double
    SphericalConstant<dim>::
    boundary_composition (const types::boundary_id boundary_indicator,
                          const Point<dim> &/*position*/,
                          const unsigned int compositional_field) const
    {
      AssertThrow (boundary_indicator == outer_boundary_indicator || boundary_indicator == inner_boundary_indicator,
                   ExcMessage ("Unknown boundary indicator for geometry model. "
                               "The given boundary should be ``top'' or ``bottom''."));

      // In case not all fields are fixed on the boundary,
      // figure out the right index of the given field.
      unsigned int field_id = compositional_field;
      if (this->get_boundary_composition_manager().boundaries_with_fixed_subset_of_fields_exist())
        {
          const std::set<unsigned int> fixed_fields = this->get_boundary_composition_manager().get_fixed_compositional_fields_for_plugin_on_boundary("spherical constant", boundary_indicator);
          Assert (fixed_fields.find(compositional_field) != fixed_fields.end(),
                  ExcMessage ("Boundary composition was requested for field " +
                              Utilities::int_to_string(compositional_field) +
                              " on boundary " +
                              Utilities::int_to_string(boundary_indicator) +
                              " but this field is not prescribed by the `spherical constant` plugin. "));
          field_id = std::distance(fixed_fields.begin(), fixed_fields.find(compositional_field));
        }

      if (boundary_indicator == outer_boundary_indicator)
        return outer_composition[field_id];
      else if (boundary_indicator == inner_boundary_indicator)
        return inner_composition[field_id];

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
                             "one entry or as many entries as there are compositional fields "
                             "prescribed by the plugin on the boundary. Units: none.");
          prm.declare_entry ("Inner composition", "1.",
                             Patterns::List(Patterns::Double ()),
                             "A comma separated list of composition boundary values "
                             "at the bottom boundary (at minimal radius). This list must have "
                             "one entry or as many entries as there are compositional fields "
                             "prescribed by the plugin on the boundary. Units: none.");
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
          // Set the boundary id of the inner boundary and get those compositional
          // fields that are fixed on it. There is no inner boundary in a full sphere.
          // Assume that if fields are only fixed on one boundary, no values are set
          // for the other boundary. Therefore the default one input value will be used
          // for possibly_extend_from_1_to_N, which will return a vector of size one,
          // without comparing to N = 0.
          unsigned int n_inner_fixed_fields = 0;
          if (Plugins::plugin_type_matches<const GeometryModel::Sphere<dim>>(this->get_geometry_model()))
            inner_boundary_indicator = numbers::invalid_unsigned_int;
          else
            {
              inner_boundary_indicator = this->get_geometry_model().translate_symbolic_boundary_name_to_id("bottom");
              n_inner_fixed_fields = (this->get_boundary_composition_manager().get_fixed_compositional_fields_for_plugin_on_boundary("spherical constant", inner_boundary_indicator)).size();
            }

          // Get the boundary id of the outer boundary, which exists in all spherical geometries.
          outer_boundary_indicator = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");

          // Get the compositional fields that are fixed on outer boundary.
          const unsigned int n_outer_fixed_fields = (this->get_boundary_composition_manager().get_fixed_compositional_fields_for_plugin_on_boundary("spherical constant", outer_boundary_indicator)).size();

          Assert (n_inner_fixed_fields <= n_compositional_fields && n_outer_fixed_fields <= n_compositional_fields,
                  ExcMessage ("The number of fixed compositional fields on the inner and/or outer boundary is higher than the total number of fields."));

          // Read in the composition values (either one value or as many as there are fixed fields).
          inner_composition = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Inner composition"))),
                                                                      n_inner_fixed_fields,
                                                                      "Inner boundary composition values");

          outer_composition = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Outer composition"))),
                                                                      n_outer_fixed_fields,
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
