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


#include <aspect/boundary_composition/two_merged_boxes.h>
#include <aspect/geometry_model/two_merged_boxes.h>

#include <utility>
#include <limits>


namespace aspect
{
  namespace BoundaryComposition
  {

    template <int dim>
    double
    TwoMergedBoxes<dim>::
    boundary_composition (const types::boundary_id boundary_indicator,
                          const Point<dim> &/*position*/,
                          const unsigned int compositional_field) const
    {
      // verify that the geometry is in fact a box since only
      // for this geometry do we know for sure what boundary indicators it
      // uses and what they mean
      Assert (dynamic_cast<const GeometryModel::TwoMergedBoxes<dim>*>(&this->get_geometry_model())
              != 0,
              ExcMessage ("This boundary model is only useful if the geometry is "
                          "in fact a box with additional lithosphere boundary indicators."));

      Assert (boundary_indicator<2*dim+2*(dim-1), ExcMessage ("The given boundary indicator needs to be less than 2*dimension+2*(dim-1)."));

      return composition_values[boundary_indicator][compositional_field];
    }

    template <int dim>
    void
    TwoMergedBoxes<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary composition model");
      {
        prm.enter_subsection("Box with lithosphere boundary indicators");
        {
          prm.declare_entry ("Left composition", "",
                             Patterns::List(Patterns::Double ()),
                             "A comma separated list of composition boundary values "
                             "at the left boundary (at minimal x-value). This list must have as many "
                             "entries as there are compositional fields. Units: none.");
          prm.declare_entry ("Right composition", "",
                             Patterns::List(Patterns::Double ()),
                             "A comma separated list of composition boundary values "
                             "at the right boundary (at maximal x-value). This list must have as many "
                             "entries as there are compositional fields. Units: none.");
          prm.declare_entry ("Left composition lithosphere", "",
                             Patterns::List(Patterns::Double ()),
                             "A comma separated list of composition boundary values "
                             "at the left boundary (at minimal x-value). This list must have as many "
                             "entries as there are compositional fields. Units: none.");
          prm.declare_entry ("Right composition lithosphere", "",
                             Patterns::List(Patterns::Double ()),
                             "A comma separated list of composition boundary values "
                             "at the right boundary (at maximal x-value). This list must have as many "
                             "entries as there are compositional fields. Units: none.");
          prm.declare_entry ("Bottom composition", "",
                             Patterns::List(Patterns::Double ()),
                             "A comma separated list of composition boundary values "
                             "at the bottom boundary (at minimal y-value in 2d, or minimal "
                             "z-value in 3d). This list must have as many "
                             "entries as there are compositional fields. Units: none.");
          prm.declare_entry ("Top composition", "",
                             Patterns::List(Patterns::Double ()),
                             "A comma separated list of composition boundary values "
                             "at the top boundary (at maximal y-value in 2d, or maximal "
                             "z-value in 3d). This list must have as many "
                             "entries as there are compositional fields. Units: none.");
          if (dim==3)
            {
              prm.declare_entry ("Front composition", "",
                                 Patterns::List(Patterns::Double ()),
                                 "A comma separated list of composition boundary values "
                                 "at the front boundary (at minimal y-value). This list must have as many "
                                 "entries as there are compositional fields. Units: none.");
              prm.declare_entry ("Back composition", "",
                                 Patterns::List(Patterns::Double ()),
                                 "A comma separated list of composition boundary values "
                                 "at the back boundary (at maximal y-value). This list must have as many "
                                 "entries as there are compositional fields. Units: none.");
              prm.declare_entry ("Front composition lithosphere", "",
                                 Patterns::List(Patterns::Double ()),
                                 "A comma separated list of composition boundary values "
                                 "at the front lithosphere boundary (at momimal y-value). This list must have as many "
                                 "entries as there are compositional fields. Units: none.");
              prm.declare_entry ("Back composition lithosphere", "",
                                 Patterns::List(Patterns::Double ()),
                                 "A comma separated list of composition boundary values "
                                 "at the back lithosphere boundary (at maximal y-value). This list must have as many "
                                 "entries as there are compositional fields. Units: none.");
            }
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    TwoMergedBoxes<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary composition model");
      {
        prm.enter_subsection("Box with lithosphere boundary indicators");
        {
          switch (dim)
            {
              case 2:
                composition_values[0] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Left composition")));
                composition_values[1] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Right composition")));
                composition_values[2] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Bottom composition")));
                composition_values[3] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Top composition")));
                composition_values[4] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Left composition lithosphere")));
                composition_values[5] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Right composition lithosphere")));
                break;

              case 3:
                composition_values[0] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Left composition")));
                composition_values[1] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Right composition")));
                composition_values[2] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Front composition")));
                composition_values[3] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Back composition")));
                composition_values[4] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Bottom composition")));
                composition_values[5] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Top composition")));
                composition_values[6] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Left composition lithosphere")));
                composition_values[7] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Right composition lithosphere")));
                composition_values[8] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Front composition lithosphere")));
                composition_values[9] = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Back composition lithosphere")));
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
    TwoMergedBoxes<dim>::initialize()
    {
      // verify that each of the lists for boundary values
      // has the requisite number of elements
      for (unsigned int f=0; f<2*dim+2*(dim-1); ++f)
        AssertThrow (composition_values[f].size() == this->n_compositional_fields(),
                     ExcMessage (std::string("The specification of boundary composition values for the 'box with "
                                             "lithosphere boundary indicators' model "
                                             "requires as many values for each boundary indicatof the box as there "
                                             "are compositional "
                                             "fields. However, for boundary indicator ")
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
    ASPECT_REGISTER_BOUNDARY_COMPOSITION_MODEL(TwoMergedBoxes,
                                               "box with lithosphere boundary indicators",
                                               "A model in which the composition is chosen constant on "
                                               "all the sides of a box. Additional boundary indicators "
                                               "are added to the lithospheric parts of the vertical boundaries. "
                                               "This model is to be used with the 'Two Merged Boxes' Geometry Model.")
  }
}
