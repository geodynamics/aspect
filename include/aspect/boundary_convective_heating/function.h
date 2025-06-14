/*
  Copyright (C) 2025 - by the authors of the ASPECT code.

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


#ifndef _aspect_boundary_convective_heating_function_h
#define _aspect_boundary_convective_heating_function_h

#include <aspect/boundary_convective_heating/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace BoundaryConvectiveHeating
  {
    /**
     * A class that implements boundary convective heating (Robin) boundary
     * conditions based on a functional description provided in the input file.
     *
     * @ingroup BoundaryConvectiveHeating
     */
    template <int dim>
    class Function : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Compute the heat transfer coefficients for a list of evaluation points.
         *
         * @param boundary_indicator The boundary indicator of the part of the
         * boundary of the domain on which the evaluation points are located
         * and where we are requesting the heat transfer coefficients.
         * @param material_model_inputs The material property inputs.
         * @param material_model_outputs The material property outputs.
         *
         * @return A vector of heat transfer coefficients at the evaluation points.
         */
        std::vector<double>
        heat_transfer_coefficient (const types::boundary_id boundary_indicator,
                                   const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                                   const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs) const override;

        /**
         * A function that is called at the beginning of each time step to
         * indicate what the model time is for which the boundary values will
         * next be evaluated. For the current class, the function passes to
         * the parsed function what the current time is.
         */
        void update () override;

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
         * A function object representing the boundary convective heating.
         */
        Functions::ParsedFunction<dim> boundary_convective_heating_function;

        /**
         * The coordinate representation to evaluate the function. Possible
         * choices are depth, cartesian and spherical.
         */
        Utilities::Coordinates::CoordinateSystem coordinate_system;
    };
  }
}


#endif
