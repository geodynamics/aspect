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


#include <aspect/gravity_model/vertical.h>

#include <deal.II/base/tensor.h>

namespace aspect
{
  namespace GravityModel
  {
    template <int dim>
    Tensor<1,dim>
    Vertical<dim>::gravity_vector (const Point<dim> &) const
    {
      Tensor<1,dim> g;
      g[dim-1] = -1;
      return g;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace GravityModel
  {
    ASPECT_REGISTER_GRAVITY_MODEL(Vertical,
                                  "vertical",
                                  "A gravity model in which the gravity direction is vertically downward "
                                  "and at a constant magnitude equal to one.");
  }
}
