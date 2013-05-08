/*
  Copyright (C) 2013 by the authors of the ASPECT code.

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
/*  $Id: box.cc 1545 2013-01-21 16:14:41Z heister $  */


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
    composition (const GeometryModel::Interface<dim> &geometry_model,
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
      return composition_[boundary_indicator];
    }


    template <int dim>
    double
    Box<dim>::
    minimal_composition (const std::set<types::boundary_id>& fixed_boundary_ids) const
    {
      if (fixed_boundary_ids.empty())
    	return *std::min_element(composition_, composition_+2*dim);
      else
      {
        double min = maximal_composition(fixed_boundary_ids);
        for (typename std::set<types::boundary_id>::const_iterator
             p = fixed_boundary_ids.begin();
             p != fixed_boundary_ids.end(); ++p)
          if (p != fixed_boundary_ids.end())
            min = std::min(min,composition_[*p]);
        return min;
      }
    }



    template <int dim>
    double
    Box<dim>::
    maximal_composition (const std::set<types::boundary_id>& fixed_boundary_ids) const
    {
      if (fixed_boundary_ids.empty())
    	return *std::max_element(composition_, composition_+2*dim);
      else
      {
    	double max = std::numeric_limits<double>::min();
        for (typename std::set<types::boundary_id>::const_iterator
             p = fixed_boundary_ids.begin();
             p != fixed_boundary_ids.end(); ++p)
          if (p != fixed_boundary_ids.end())
            max = std::max(max,composition_[*p]);
        return max;
      }
    }

    template <int dim>
    void
    Box<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary composition model");
      {
        prm.enter_subsection("Box");
        {
          prm.declare_entry ("Left composition", "0",
                             Patterns::Double (),
                             "Composition at the left boundary (at minimal x-value). Units: K.");
          prm.declare_entry ("Right composition", "0",
                             Patterns::Double (),
                             "Composition at the right boundary (at maximal x-value). Units: K.");
          prm.declare_entry ("Bottom composition", "1",
                             Patterns::Double (),
                             "Composition at the bottom boundary (at minimal z-value). Units: K.");
          prm.declare_entry ("Top composition", "0",
                             Patterns::Double (),
                             "Composition at the top boundary (at maximal x-value). Units: K.");
          if (dim==3)
            {
              prm.declare_entry ("Front composition", "0",
                                 Patterns::Double (),
                                 "Composition at the front boundary (at minimal y-value). Units: K.");
              prm.declare_entry ("Back composition", "0",
                                 Patterns::Double (),
                                 "Composition at the back boundary (at maximal y-value). Units: K.");
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
                composition_[0] = prm.get_double ("Left composition");
                composition_[1] = prm.get_double ("Right composition");
                composition_[2] = prm.get_double ("Bottom composition");
                composition_[3] = prm.get_double ("Top composition");
                break;

              case 3:
                composition_[0] = prm.get_double ("Left composition");
                composition_[1] = prm.get_double ("Right composition");
                composition_[2] = prm.get_double ("Front composition");
                composition_[3] = prm.get_double ("Back composition");
                composition_[4] = prm.get_double ("Bottom composition");
                composition_[5] = prm.get_double ("Top composition");
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
  namespace BoundaryComposition
  {
    ASPECT_REGISTER_BOUNDARY_COMPOSITION_MODEL(Box,
                                               "box",
                                               "A model in which the composition is chosen constant on "
                                               "all the sides of a box.")
  }
}
