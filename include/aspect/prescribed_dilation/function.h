/*
  Copyright (C) 2014 - 2024 by the authors of the ASPECT code.

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


#ifndef _aspect_prescribed_dilation_function_h
#define _aspect_prescribed_dilation_function_h

#include <aspect/prescribed_dilation/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace PrescribedDilation
  {
    /**
     * A class that implements a dilation model based on a functional
     * description provided in the input file.
     *
     * @ingroup PrescribedDilation
     */
    template <int dim>
    class Function : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        Function();

        /**
         * Compute dilation based on the function and store the result in the PrescribedDilationOutputs structure
         */
        void
        evaluate(const aspect::MaterialModel::MaterialModelInputs<dim> &in, PrescribedDilationOutputs<dim> &out) const override;

        /**
         * A function that is called at the beginning of each time step to
         * allow the model to do whatever necessary. In this case the time of
         * the function object is updated.
         */
        void
        update() override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static void
        declare_parameters(ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters(ParameterHandler &prm) override;

      private:
        /**
         * A function object representing the dilation.
         */
        Functions::ParsedFunction<dim> function;

        /**
         * The coordinate representation to evaluate the function. Possible
         * choices are depth, cartesian and spherical.
         */
        Utilities::Coordinates::CoordinateSystem coordinate_system;
    };
  }
}


#endif
