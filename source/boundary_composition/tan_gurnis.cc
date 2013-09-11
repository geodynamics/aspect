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
/*  $Id: tan_gurnis.cc 1538 2013-01-06 03:12:23Z bangerth $  */


#include <aspect/boundary_composition/tan_gurnis.h>
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
    TanGurnis<dim>::
    composition (const GeometryModel::Interface<dim> &geometry_model,
                 const unsigned int                   boundary_indicator,
                 const Point<dim>                    &location,
                 const unsigned int                   compositional_field) const
    {
      // verify that the geometry is in fact a box since only
      // for this geometry do we know for sure what boundary indicators it
      // uses and what they mean
      Assert (dynamic_cast<const GeometryModel::Box<dim>*>(&geometry_model)
              != 0,
              ExcMessage ("This boundary model is only implemented if the geometry is "
                          "in fact a box."));

      double wavenumber=1;
      return sin(numbers::PI*location(dim-1))*cos(numbers::PI*wavenumber*location(0));
    }


    template <int dim>
    double
    TanGurnis<dim>::
    minimal_composition (const std::set<types::boundary_id>& fixed_boundary_ids) const
    {
      return 0;
    }



    template <int dim>
    double
    TanGurnis<dim>::
    maximal_composition (const std::set<types::boundary_id>& fixed_boundary_ids) const
    {
      return 1;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryComposition
  {
    ASPECT_REGISTER_BOUNDARY_COMPOSITION_MODEL(TanGurnis,
                                               "Tan Gurnis",
                                               "A model for the Tan/Gurnis benchmark.")
  }
}
