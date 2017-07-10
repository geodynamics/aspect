/*
  Copyright (C) 2014 - 2016 by the authors of the ASPECT code.

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


#ifndef _aspect_gravity_model_radial_h
#define _aspect_gravity_model_radial_h

#include <aspect/simulator_access.h>
#include <aspect/gravity_model/interface.h>

namespace aspect
{
  namespace GravityModel
  {
    using namespace dealii;

    /**
     * A class that describes gravity as a radial vector of constant
     * magnitude. The magnitude's value is read from the input file.
     *
     * @ingroup GravityModels
     */
    template <int dim>
    class RadialConstant : public Interface<dim>
    {
      public:
        /**
         * Return the gravity vector as a function of position.
         */
        virtual Tensor<1,dim> gravity_vector (const Point<dim> &position) const;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        /**
         * Magnitude of the gravity vector.
         */
        double magnitude;
    };


    /**
     * A class that describes gravity as a radial vector with a magnitude that
     * is physically correct for an earth in which there are different
     * densities for the earth core and mantle. Specifically, at the core-
     * mantle boundary, gravity is assumed to be equal to 10.7 m/s^2 and it is
     * 9.8 at the earth surface; in between, it follows the behavior one would
     * expect for a mantle of constant density.
     *
     * This is the model used and discussed in the step-32 tutorial program of
     * deal.II.
     *
     * @ingroup GravityModels
     */
    template <int dim>
    class RadialEarthLike : public Interface<dim>
    {
      public:
        /**
         * Return the gravity vector as a function of position.
         */
        virtual Tensor<1,dim> gravity_vector (const Point<dim> &position) const;
    };


    /**
     * A class that describes gravity as a radial vector of linearly
     * decreasing magnitude with depth.  Meant for use in the Sphere geometry
     * model, where you expect that kind of field if one assumed a constant
     * density of Earth.
     *
     * @ingroup GravityModels
     */
    template <int dim>
    class RadialLinear : public Interface<dim>, public virtual SimulatorAccess<dim>
    {
      public:
        /**
         * Return the gravity vector as a function of position.
         */
        virtual Tensor<1,dim> gravity_vector (const Point<dim> &position) const;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        /**
         * Magnitude of the gravity vector at the surface, m/s^2
         */
        double magnitude_at_surface;

    };
  }
}

#endif
