// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2022-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

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
