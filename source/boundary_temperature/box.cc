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


#include <aspect/boundary_temperature/box.h>
#include <aspect/geometry_model/box.h>

#include <utility>
#include <limits>


namespace aspect
{
  namespace BoundaryTemperature
  {
// ------------------------------ Box -------------------

    template <int dim>
    double
    Box<dim>::
    temperature (const GeometryModel::Interface<dim> &geometry_model,
                 const unsigned int                   boundary_indicator,
                 const Point<dim>                    &location) const
    {
      // verify that the geometry is in fact a box since only
      // for this geometry do we know for sure what boundary indicators it
      // uses and what they mean
      Assert (dynamic_cast<const GeometryModel::Box<dim>*>(&geometry_model)
              != 0,
              ExcMessage ("This boundary model is only implemented if the geometry is "
                          "in fact a box."));

      Assert (boundary_indicator<2*dim, ExcMessage ("Unknown boundary indicator."));
      return temperature_[boundary_indicator];
    }


    template <int dim>
    double
    Box<dim>::
    minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      if (fixed_boundary_ids.empty())
        return *std::min_element(temperature_, temperature_+2*dim);
      else
        {
          double min = maximal_temperature(fixed_boundary_ids);
          for (typename std::set<types::boundary_id>::const_iterator
               p = fixed_boundary_ids.begin();
               p != fixed_boundary_ids.end(); ++p)
            if (p != fixed_boundary_ids.end())
              min = std::min(min,temperature_[*p]);
          return min;
        }
    }



    template <int dim>
    double
    Box<dim>::
    maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      if (fixed_boundary_ids.empty())
        return *std::max_element(temperature_, temperature_+2*dim);
      else
        {
          double max = -std::numeric_limits<double>::max();
          for (typename std::set<types::boundary_id>::const_iterator
               p = fixed_boundary_ids.begin();
               p != fixed_boundary_ids.end(); ++p)
            if (p != fixed_boundary_ids.end())
              max = std::max(max,temperature_[*p]);
          return max;
        }
    }

    template <int dim>
    void
    Box<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary temperature model");
      {
        prm.enter_subsection("Box");
        {
          prm.declare_entry ("Left temperature", "1",
                             Patterns::Double (),
                             "Temperature at the left boundary (at minimal x-value). Units: K.");
          prm.declare_entry ("Right temperature", "0",
                             Patterns::Double (),
                             "Temperature at the right boundary (at maximal x-value). Units: K.");
          prm.declare_entry ("Bottom temperature", "0",
                             Patterns::Double (),
                             "Temperature at the bottom boundary (at minimal z-value). Units: K.");
          prm.declare_entry ("Top temperature", "0",
                             Patterns::Double (),
                             "Temperature at the top boundary (at maximal x-value). Units: K.");
          if (dim==3)
            {
              prm.declare_entry ("Front temperature", "0",
                                 Patterns::Double (),
                                 "Temperature at the front boundary (at minimal y-value). Units: K.");
              prm.declare_entry ("Back temperature", "0",
                                 Patterns::Double (),
                                 "Temperature at the back boundary (at maximal y-value). Units: K.");
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
      prm.enter_subsection("Boundary temperature model");
      {
        prm.enter_subsection("Box");
        {
          switch (dim)
            {
              case 2:
                temperature_[0] = prm.get_double ("Left temperature");
                temperature_[1] = prm.get_double ("Right temperature");
                temperature_[2] = prm.get_double ("Bottom temperature");
                temperature_[3] = prm.get_double ("Top temperature");
                break;

              case 3:
                temperature_[0] = prm.get_double ("Left temperature");
                temperature_[1] = prm.get_double ("Right temperature");
                temperature_[2] = prm.get_double ("Front temperature");
                temperature_[3] = prm.get_double ("Back temperature");
                temperature_[4] = prm.get_double ("Bottom temperature");
                temperature_[5] = prm.get_double ("Top temperature");
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
    ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL(Box,
                                               "box",
                                               "A model in which the temperature is chosen constant on "
                                               "all the sides of a box.")
  }
}
