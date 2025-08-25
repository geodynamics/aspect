/*
  Copyright (C) 2011 - 2025 by the authors of the ASPECT code.

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

#ifndef _aspect_prescribed_solution_velocity_function_h
#define _aspect_prescribed_solution_velocity_function_h

#include <aspect/prescribed_solution/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace PrescribedSolution
  {

    /**
     * A class that implements prescribed fields conditions based on a
     * functional description provided in the input file.
     *
     * @ingroup PrescribedSolution
     */
    template <int dim>
    class VelocityFunction : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:

        /**
         * Constructor.
         */
        VelocityFunction ();

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
         * Decide and assign cell-wise constraints for velocity DoFs.
         * This function inspects component indices of constraints at local
         * DoFs in a cell. An indicator function is evaluated at every position
         * to decide whether new constraints need to be added. If so, values of
         * new constraints are computed and stored in the solution.
         */
        void constrain_solution (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                 const std::vector<Point<dim>> &positions,
                                 const std::vector<unsigned int> component_indices,
                                 std::vector<bool> &should_be_constrained,
                                 std::vector<double> &solution) override;

      private:
        /**
         * A function object representing the indicator function for prescribed solution
         */
        Functions::ParsedFunction<dim> prescribed_velocity_indicator_function;

        /**
         * A function object representing the components of the velocity.
         */
        Functions::ParsedFunction<dim> prescribed_velocity_function;

        /**
         * The coordinate representation to evaluate the function. Possible
         * choices are depth, cartesian and spherical.
         */
        Utilities::Coordinates::CoordinateSystem coordinate_system;

        /**
         * Whether to specify velocity in x, y, z components, or
         * r, phi, theta components.
         */
        bool use_spherical_unit_vectors;
    };
  }
}


#endif
