/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
/*  $Id$  */


#include <aspect/boundary_temperature/tan_gurnis.h>
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
    TanGurnis<dim>::
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

      double wavenumber=1;
      return sin(numbers::PI*location(dim-1))*cos(numbers::PI*wavenumber*location(0));
    }


    template <int dim>
    double
    TanGurnis<dim>::
    minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      return 0;
    }



    template <int dim>
    double
    TanGurnis<dim>::
    maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      return 1;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryTemperature
  {
    ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL(TanGurnis,
                                               "Tan Gurnis",
                                               "A model for the Tan/Gurnis benchmark.")
  }
}
