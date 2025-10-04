/*
  Copyright (C) 2014 - 2019 by the authors of the ASPECT code.

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


#ifndef _aspect_gravity_model_radial_earth_like_h
#define _aspect_gravity_model_radial_earth_like_h

#include <aspect/simulator_access.h>
#include <aspect/gravity_model/interface.h>

namespace aspect
{
  namespace GravityModel
  {
    /**
     * This model has been removed due to its misleading name. The available
     * AsciiData gravity model (using default parameters) is much more
     * earth-like, since it uses the gravity profile used in the construction
     * of the Preliminary Reference Earth Model (PREM, Dziewonski and Anderson,
     * 1981).
     *
     * This is the model used and discussed in the step-32 tutorial program of
     * deal.II.
     *
     * @ingroup GravityModels
     */
    template <int dim>
    class RadialEarthLike : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Initialization.
         */
        void initialize() override;

        /**
         * Return the gravity vector as a function of position.
         */
        Tensor<1,dim> gravity_vector (const Point<dim> &position) const override;
    };
  }
}

#endif
