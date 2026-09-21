// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2014-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_boundary_composition_function_h
#define _aspect_boundary_composition_function_h

#include <aspect/boundary_composition/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace BoundaryComposition
  {
    /**
     * A class that implements boundary composition based on a functional
     * description provided in the input file.
     *
     * @ingroup BoundaryCompositions
     */
    template <int dim>
    class Function : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the boundary composition as a function of position and time.
         *
         * @copydoc aspect::BoundaryComposition::Interface::boundary_composition()
         */
        double boundary_composition (const types::boundary_id boundary_indicator,
                                     const Point<dim> &position,
                                     const unsigned int compositional_field) const override;
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
        void declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * A function object representing the compositional fields.
         */
        std::unique_ptr<Functions::ParsedFunction<dim>> function;

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
