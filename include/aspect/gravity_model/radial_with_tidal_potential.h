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


#ifndef _aspect_gravity_model_radial_with_tidal_potential_h
#define _aspect_gravity_model_radial_with_tidal_potential_h

#include <aspect/simulator_access.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/gravity_model/radial.h>

namespace aspect
{
  namespace GravityModel
  {
    /**
     * A class that describes gravity as a radial vector of linearly
     * changing magnitude, which is modified by a tidal potential from flattening and non-synchronous rotation.
     *
     * The equation implemented in this gravity model is from Tobie et al. (2025) (https://doi.org/10.1007/s11214-025-01136-y),
     * which is defined as:
     * g = -magnitude - gradient (-(tidal potential)).
     * Tidal potential is positive because the formula follows conventions from geodesy research, where potential is taken as positive.
     * (tidal potential) = (3 G M_p) / (2 a_s^3) * r^2 * (Tstar + T0)
     * Tstar = 1/6 *(1-3*cos(theta)^2) and T0=1/2sin(theta)^2*cos(2*lambda + 2*b*t)
     * where G = gravitational constant, M_p = mass of the perturbing body, a_s = semimajor axis of the orbit, b = angular rate of non-synchronous rotation.
     * b = 2 * pi / P where P is period of NSR.
     * r, theta and lambda are radial distance, polar angle and azimuthal angle, respectively.
     *
     * @ingroup GravityModels
     */
    template <int dim>
    class RadialWithTidalPotential : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the gravity vector as a function of position.
         */
        Tensor<1,dim> gravity_vector (const Point<dim> &position) const override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * Mass of the perturbing body
         */
        double M_p;

        /**
         * Semimajor axis of the orbit that causes the tidal perturbation
         */
        double a_s;

        /**
         * Period of the non-synchronous rotation in year or second
         */
        double P;
    };
  }
}

#endif
