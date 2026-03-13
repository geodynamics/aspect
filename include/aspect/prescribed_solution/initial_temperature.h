/*
  Copyright (C) 2011 - 2026 by the authors of the ASPECT code.

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
  along with ASPECT; see the file LICENSE. If not see
  <http://www.gnu.org/licenses/>.
*/

#ifndef _aspect_prescribed_solution_initial_temperature_h
#define _aspect_prescribed_solution_initial_temperature_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/prescribed_solution/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace PrescribedSolution
  {
    /**
     * A class that prescribes temperature in a selected region using
     * the value returned by the active initial temperature model.
     *
     * The region is selected through an indicator function. Where the
     * indicator is greater than 0.5, the temperature is constrained
     * to the initial temperature value evaluated at that position.
     *
     * @ingroup PrescribedSolution
     */
    template <int dim>
    class InitialTemperature
      : public Interface<dim>,
        public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        InitialTemperature ();

        /**
        * Store a shared pointer to the initial temperature manager so the
        * plugin can safely access it after initialization.
        */
        void initialize () override;

        /**
         * Update the current time in the indicator function.
         */
        void update () override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static void declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void parse_parameters (ParameterHandler &prm) override;

        /**
         * Decide and assign cell-wise constraints for temperature DoFs.
         */
        void constrain_solution (
          const typename DoFHandler<dim>::active_cell_iterator &cell,
          const std::vector<Point<dim>> &positions,
          const std::vector<unsigned int> &component_indices,
          std::vector<bool> &should_be_constrained,
          std::vector<double> &solution) override;

      private:
        /**
         * Indicator function for selecting where the temperature should
         * be prescribed.
         */
        Functions::ParsedFunction<dim> prescribed_temperature_indicator_function;

        /**
         * The coordinate representation used to evaluate the indicator
         * function. Possible choices are cartesian, spherical, and depth.
         */
        Utilities::Coordinates::CoordinateSystem coordinate_system;

        /**
         * Shared pointer to the initial temperature manager. We keep this
         * alive because the simulator may release its own pointer after
         * initialization.
         */
        std::shared_ptr<const aspect::InitialTemperature::Manager<dim>> initial_temperature_manager;
    };
  }
}

#endif
