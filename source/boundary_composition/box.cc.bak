/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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


#include <aspect/boundary_composition/box.h>
#include <aspect/geometry_model/box.h>

#include <utility>
#include <limits>


namespace aspect
{
  namespace BoundaryComposition
  {
// ------------------------------ Box -------------------

    template <int dim>
    double
    Box<dim>::
    boundary_composition (const types::boundary_id boundary_indicator,
                          const Point<dim> &/*position*/,
                          const unsigned int compositional_field) const
    {
      Assert (boundary_indicator<2*dim, ExcMessage ("The given boundary indicator needs to be less than 2*dimension.."));
      return composition_values[boundary_indicator][compositional_field];
    }

    template <int dim>
    void
    Box<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary composition model");
      {
        prm.enter_subsection("Box");
        {
          prm.declare_entry ("Left composition", "",
                             Patterns::List(Patterns::Double ()),
                             "A comma separated list of composition boundary values "
                             "at the left boundary (at minimal $x$-value). This list must have as many "
                             "entries as there are compositional fields. Units: none.");
          prm.declare_entry ("Right composition", "",
                             Patterns::List(Patterns::Double ()),
                             "A comma separated list of composition boundary values "
                             "at the right boundary (at maximal $x$-value). This list must have as many "
                             "entries as there are compositional fields. Units: none.");
          prm.declare_entry ("Bottom composition", "",
                             Patterns::List(Patterns::Double ()),
                             "A comma separated list of composition boundary values "
                             "at the bottom boundary (at minimal $y$-value in 2d, or minimal "
                             "$z$-value in 3d). This list must have as many "
                             "entries as there are compositional fields. Units: none.");
          prm.declare_entry ("Top composition", "",
                             Patterns::List(Patterns::Double ()),
                             "A comma separated list of composition boundary values "
                             "at the top boundary (at maximal $y$-value in 2d, or maximal "
                             "$z$-value in 3d). This list must have as many "
                             "entries as there are compositional fields. Units: none.");
          if (dim==3)
            {
              prm.declare_entry ("Front composition", "",
                                 Patterns::List(Patterns::Double ()),
                                 "A comma separated list of composition boundary values "
                                 "at the front boundary (at minimum $y$-value). This list must have as many "
                                 "entries as there are compositional fields. Units: none.");
              prm.declare_entry ("Back composition", "",
                                 Patterns::List(Patterns::Double ()),
                                 "A comma separated list of composition boundary values "
                                 "at the back boundary (at maximum $y$-value). This list must have as many "
                                 "entries as there are compositional fields. Units: none.");
            }
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    Box<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary composition model");
      {
        prm.enter_subsection("Box");
        {
          switch (dim)
            {
              case 2:
                composition_values[0] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Left composition")));
                composition_values[1] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Right composition")));
                composition_values[2] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Bottom composition")));
                composition_values[3] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Top composition")));
                break;

              case 3:
                composition_values[0] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Left composition")));
                composition_values[1] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Right composition")));
                composition_values[2] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Front composition")));
                composition_values[3] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Back composition")));
                composition_values[4] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Bottom composition")));
                composition_values[5] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Top composition")));
                break;

              default:
                Assert (false, ExcNotImplemented());
            }
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }



    template <int dim>
    void
    Box<dim>::initialize()
    {
      // verify that the geometry is a box since only for this geometry
      // do we know for sure what boundary indicators it uses and what they mean
      AssertThrow (Plugins::plugin_type_matches<const GeometryModel::Box<dim>>(this->get_geometry_model()),
                   ExcMessage ("This boundary model is only implemented if the geometry is "
                               "a box."));

      // Verify that each of the lists for boundary values
      // has the requisite number of elements if it is in the set
      // of prescribed boundary indicators.
      for (unsigned int f=0; f<2*dim; ++f)
        if (this->get_boundary_composition_manager().get_fixed_composition_boundary_indicators().count(f) != 0)
          AssertThrow (composition_values[f].size() == this->n_compositional_fields(),
                       ExcMessage (std::string("The specification of boundary composition values for the `box' model "
                                               "requires as many values on each face of the box as there are compositional "
                                               "fields. However, for face ")
                                   +
                                   Utilities::int_to_string(f)
                                   +
                                   ", the input file specifies "
                                   +
                                   Utilities::int_to_string(composition_values[f].size())
                                   +
                                   " values even though there are "
                                   +
                                   Utilities::int_to_string(this->n_compositional_fields())
                                   +
                                   " compositional fields."));
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryComposition
  {
    ASPECT_REGISTER_BOUNDARY_COMPOSITION_MODEL(Box,
                                               "box",
                                               "A model in which the composition is chosen constant on "
                                               "the sides of a box which are selected by the parameters "
                                               "Left/Right/Top/Bottom/Front/Back composition")
  }
}
