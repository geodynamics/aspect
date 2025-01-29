/*
  Copyright (C) 2022 by the authors of the ASPECT code.

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

#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/gravity_model/interface.h>

namespace aspect
{
  template <int dim>
  class MyGravity :
    public aspect::GravityModel::Interface<dim>
  {
    public:
      virtual Tensor<1,dim> gravity_vector (const Point<dim> &position) const
      {
        Tensor<1,dim> ret;
        ret[0] = position[1];
        ret[1] = 42.0;
        return ret;
      }
  };


// explicit instantiation
  ASPECT_REGISTER_GRAVITY_MODEL(MyGravity, "my gravity", "no description")
}
