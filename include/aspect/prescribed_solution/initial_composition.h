/*
 Copyright (C) 2015 - 2022 by the authors of the ASPECT code.

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

#ifndef _aspect_prescribed_solution_initial_composition_h
#define _aspect_prescribed_solution_initial_composition_h

#include <aspect/initial_composition/interface.h>
#include <aspect/prescribed_solution/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace PrescribedSolution
  {
    /**
     * Prescribe compositional fields using the initial composition model.
     * The region where the constraint applies is determined by an
     * indicator function.
     */
    template <int dim>
    class InitialComposition
      : public Interface<dim>,
        public SimulatorAccess<dim>
    {
      public:

        InitialComposition ();

        /**
         * Store shared pointer to initial composition manager.
         */
        void initialize () override;

        /**
         * Update function time.
         */
        void update () override;

        /**
         * Declare parameters.
         */
        static void declare_parameters (ParameterHandler &prm);

        /**
         * Parse parameters.
         */
        void parse_parameters (ParameterHandler &prm) override;

        /**
         * Constrain compositional solution.
         */
        void constrain_solution (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                 const std::vector<Point<dim>> &positions,
                                 const std::vector<unsigned int> &component_indices,
                                 std::vector<bool> &should_be_constrained,
                                 std::vector<double> &solution) override;

      private:

        /**
         * Indicator function defining region where composition
         * is prescribed.
         */
        Functions::ParsedFunction<dim> indicator_function;

        /**
         * Coordinate system used for evaluating indicator.
         */
        Utilities::Coordinates::CoordinateSystem coordinate_system;

        /**
         * Pointer to initial composition manager.
         */
        std::shared_ptr<const aspect::InitialComposition::Manager<dim>> initial_composition_manager;
    };
  }
}

#endif
