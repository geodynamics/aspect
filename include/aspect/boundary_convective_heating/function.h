// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2025-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


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
         *
         * This variable is read from the parameter file through a parameter called 'Coordinate system'.
         */
        Utilities::Coordinates::CoordinateSystem coordinate_system;
    };
  }
}


#endif
