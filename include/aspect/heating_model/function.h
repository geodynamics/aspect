/*
  Copyright (C) 2014 by the authors of the ASPECT code.

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


#ifndef __aspect__heating_model_function_h
#define __aspect__heating_model_function_h

#include <aspect/heating_model/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace HeatingModel
  {
    using namespace dealii;

    /**
     * A class that implements a heating model based on a functional
     * description provided in the input file.
     *
     * @ingroup HeatingModels
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
         * Return the specific heating rate as calculated by the function
         * object.
         */
        virtual
        double
        specific_heating_rate (const double,
                               const double,
                               const std::vector<double> &,
                               const Point<dim> &) const;

        /**
         * A function that is called at the beginning of each time step to
         * allow the model to do whatever necessary. In this case the time of
         * the function object is updated.
         */
        virtual
        void
        update ();

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
         * A function object representing the components of the velocity.
         */
        Functions::ParsedFunction<dim> heating_model_function;
    };
  }
}


#endif
