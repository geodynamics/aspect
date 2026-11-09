// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_prescribed_solution_temperature_function_h
#define _aspect_prescribed_solution_temperature_function_h

#include <aspect/prescribed_solution/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace PrescribedSolution
  {

    /**
     * A class that prescribes the temperature in a specified region based on a
     * functional description provided in the input file.
     *
     * @ingroup PrescribedSolution
     */
    template <int dim>
    class TemperatureFunction : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:

        /**
         * Constructor.
         */
        TemperatureFunction ();

        /**
         * A function that is called at the beginning of each time step to
         * indicate what the model time is for which the solution will
         * next be prescribed. For the current class, the function passes to
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

        /**
         * Decide and assign cell-wise constraints for temperature DoFs.
         * This function inspects component indices of constraints at local
         * DoFs in a cell. An indicator function is evaluated at every position
         * to decide whether new constraints need to be added. If so, values of
         * new constraints are computed and stored in the solution.
         */
        void constrain_solution (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                 const std::vector<Point<dim>> &positions,
                                 const std::vector<unsigned int> &component_indices,
                                 std::vector<bool> &should_be_constrained,
                                 std::vector<double> &solution) override;

      private:
        /**
         * A function object representing the indicator function for prescribed temperature.
         */
        Functions::ParsedFunction<dim> prescribed_temperature_indicator_function;

        /**
         * A function object representing the prescribed temperature.
         */
        Functions::ParsedFunction<dim> prescribed_temperature_function;

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
