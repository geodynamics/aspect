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


#ifndef _aspect_gravity_model_function_h
#define _aspect_gravity_model_function_h

#include <aspect/gravity_model/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace GravityModel
  {
    using namespace dealii;

    /**
     * A class that implements gravity based on a functional description
     * provided in the input file.
     *
     * @ingroup GravityModels
     */
    template <int dim>
    class Function : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        Function ();

        /**
         * A function that is called at the beginning of each time step to
         * indicate what the model time is for which the gravity will
         * next be evaluated. For the current class, the function passes to
         * the parsed function what the current time is.
         */
        void update () override;

        /**
         * Return the gravity vector as a function of position.
         */
        Tensor<1,dim> gravity_vector (const Point<dim> &position) const override;

        /**
         * Declare the parameters this class takes through input files. The
         * default implementation of this function does not describe any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * The default implementation of this function does not read any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * A function object representing the gravity.
         */
        Functions::ParsedFunction<dim> function;
    };
  }
}


#endif
