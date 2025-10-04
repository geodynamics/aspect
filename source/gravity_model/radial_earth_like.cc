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


#include <aspect/gravity_model/radial_earth_like.h>

namespace aspect
{
  namespace GravityModel
  {
    template <int dim>
    void
    RadialEarthLike<dim>::initialize ()
    {
      AssertThrow(false,
                  ExcMessage("The 'radial earth-like' gravity model has been removed "
                             "due to its misleading name. The available AsciiData gravity "
                             "model (using default parameters) is much more earth-like, since "
                             "it uses the gravity profile used in the construction of the "
                             "Preliminary Reference Earth Model (PREM, Dziewonski and Anderson, "
                             "1981). Use the 'ascii data' model instead of 'radial earth-like'."));
    }



    template <int dim>
    Tensor<1,dim>
    RadialEarthLike<dim>::gravity_vector (const Point<dim> &/*p*/) const
    {
      // We should never get here, because of the assertion in initialize().
      AssertThrow(false,ExcInternalError());
      return Tensor<1,dim>();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace GravityModel
  {
    ASPECT_REGISTER_GRAVITY_MODEL(RadialEarthLike,
                                  "radial earth-like",
                                  "This plugin has been removed due to its misleading name. "
                                  "The included profile was hard-coded and was less earth-like "
                                  "than the `ascii data' plugin, which uses the profile "
                                  "of the Preliminary Reference Earth Model (PREM). Use `ascii data' "
                                  "instead of `radial earth-like'.")
  }
}
