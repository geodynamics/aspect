/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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


#ifndef _aspect_gravity_model_vertical_h
#define _aspect_gravity_model_vertical_h

#include <aspect/gravity_model/interface.h>

namespace aspect
{
  namespace GravityModel
  {
    using namespace dealii;

    /**
     * A class that describes gravity as a vector of constant magnitude
     * pointing vertically down.
     *
     * @ingroup GravityModels
     */
    template <int dim>
    class Vertical : public Interface<dim>
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
        double gravity_magnitude;

    };
  }
}

#endif
