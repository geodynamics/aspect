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


#include <aspect/boundary_temperature/two_merged_boxes.h>
#include <aspect/geometry_model/two_merged_boxes.h>

#include <utility>
#include <limits>


namespace aspect
{
  namespace BoundaryTemperature
  {

    template <int dim>
    double
    TwoMergedBoxes<dim>::
    boundary_temperature (const types::boundary_id boundary_indicator,
                          const Point<dim> &) const
    {
      Assert (boundary_indicator<2*dim+2*(dim-1), ExcMessage ("Given boundary indicator needs to be less than 2*dimension+2*(dimension-1)."));

      return temperature_values[boundary_indicator];
    }


    template <int dim>
    double
    TwoMergedBoxes<dim>::
    minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      if (fixed_boundary_ids.empty())
        return *std::min_element(temperature_values, temperature_values+2*dim+2*(dim-1));
      else
        {
          double min = maximal_temperature(fixed_boundary_ids);
          for (const auto p : fixed_boundary_ids)
            min = std::min(min,temperature_values[p]);
          return min;
        }
    }



    template <int dim>
    double
    TwoMergedBoxes<dim>::
    maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      if (fixed_boundary_ids.empty())
        return *std::max_element(temperature_values, temperature_values+2*dim+2*(dim-1));
      else
        {
          double max = -std::numeric_limits<double>::max();
          for (const auto p : fixed_boundary_ids)
            max = std::max(max,temperature_values[p]);
          return max;
        }
    }

    template <int dim>
    void
    TwoMergedBoxes<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary temperature model");
      {
        prm.enter_subsection("Box with lithosphere boundary indicators");
        {
          prm.declare_entry ("Left temperature", "1.",
                             Patterns::Double (),
                             "Temperature at the left boundary (at minimal $x$-value). Units: \\si{\\kelvin}.");
          prm.declare_entry ("Right temperature", "0.",
                             Patterns::Double (),
                             "Temperature at the right boundary (at maximal $x$-value). Units: \\si{\\kelvin}.");
          prm.declare_entry ("Bottom temperature", "0.",
                             Patterns::Double (),
                             "Temperature at the bottom boundary (at minimal $z$-value). Units: \\si{\\kelvin}.");
          prm.declare_entry ("Top temperature", "0.",
                             Patterns::Double (),
                             "Temperature at the top boundary (at maximal $x$-value). Units: \\si{\\kelvin}.");
          prm.declare_entry ("Left temperature lithosphere", "0.",
                             Patterns::Double (),
                             "Temperature at the additional left lithosphere boundary (specified by user in Geometry Model). Units: \\si{\\kelvin}.");
          prm.declare_entry ("Right temperature lithosphere", "0.",
                             Patterns::Double (),
                             "Temperature at the additional right lithosphere boundary (specified by user in Geometry Model). Units: \\si{\\kelvin}.");
          if (dim==3)
            {
              prm.declare_entry ("Front temperature", "0.",
                                 Patterns::Double (),
                                 "Temperature at the front boundary (at minimal $y$-value). Units: \\si{\\kelvin}.");
              prm.declare_entry ("Back temperature", "0.",
                                 Patterns::Double (),
                                 "Temperature at the back boundary (at maximal $y$-value). Units: \\si{\\kelvin}.");
              prm.declare_entry ("Front temperature lithosphere", "0.",
                                 Patterns::Double (),
                                 "Temperature at the additional front lithosphere boundary (at minimal $y$-value). Units: \\si{\\kelvin}.");
              prm.declare_entry ("Back temperature lithosphere", "0.",
                                 Patterns::Double (),
                                 "Temperature at the additional back lithosphere boundary (at maximal $y$-value). Units: \\si{\\kelvin}.");
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
      prm.enter_subsection("Boundary temperature model");
      {
        prm.enter_subsection("Box with lithosphere boundary indicators");
        {
          // verify that the geometry is a box since only for this geometry
          // do we know for sure what boundary indicators it uses and what they mean
          AssertThrow (Plugins::plugin_type_matches<const GeometryModel::TwoMergedBoxes<dim>>(this->get_geometry_model()),
                       ExcMessage ("This boundary model is only useful if the geometry is "
                                   "a box with additional boundary indicators."));

          switch (dim)
            {
              case 2:
                temperature_values[0] = prm.get_double ("Left temperature");
                temperature_values[1] = prm.get_double ("Right temperature");
                temperature_values[2] = prm.get_double ("Bottom temperature");
                temperature_values[3] = prm.get_double ("Top temperature");
                temperature_values[4] = prm.get_double ("Left temperature lithosphere");
                temperature_values[5] = prm.get_double ("Right temperature lithosphere");
                break;

              case 3:
                temperature_values[0] = prm.get_double ("Left temperature");
                temperature_values[1] = prm.get_double ("Right temperature");
                temperature_values[2] = prm.get_double ("Front temperature");
                temperature_values[3] = prm.get_double ("Back temperature");
                temperature_values[4] = prm.get_double ("Bottom temperature");
                temperature_values[5] = prm.get_double ("Top temperature");
                temperature_values[6] = prm.get_double ("Left temperature lithosphere");
                temperature_values[7] = prm.get_double ("Right temperature lithosphere");
                temperature_values[8] = prm.get_double ("Front temperature lithosphere");
                temperature_values[9] = prm.get_double ("Back temperature lithosphere");
                break;

              default:
                Assert (false, ExcNotImplemented());
            }
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
  namespace BoundaryTemperature
  {
    ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL(TwoMergedBoxes,
                                               "box with lithosphere boundary indicators",
                                               "A model in which the temperature is chosen constant on "
                                               "all the sides of a box. Additional boundary indicators "
                                               "are added to the lithospheric parts of the vertical boundaries. "
                                               "This model is to be used with the 'Two Merged Boxes' Geometry Model.")
  }
}
