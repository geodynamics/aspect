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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#ifndef _aspect_prescribed_solution_world_builder_h
#define _aspect_prescribed_solution_world_builder_h

#include <aspect/prescribed_solution/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace PrescribedSolution
  {
    /**
     * Prescribe velocity, temperature, and composition using values from the
     * WorldBuilder. WorldBuilder indicators named `temperature`, `velocity`,
     * and `composition` select the locations where the corresponding solution
     * components are constrained.
     */
    template <int dim>
    class WorldBuilder
      : public Interface<dim>,
        public SimulatorAccess<dim>
    {
      public:
        /**
         * Store a shared pointer to the WorldBuilder object.
         */
        void initialize () override;

        /**
         * Declare parameters for this plugin.
         */
        static void declare_parameters (ParameterHandler &prm);

        /**
         * Apply WorldBuilder values to solution components selected by the
         * corresponding WorldBuilder indicators.
         */
        void constrain_solution (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                 const std::vector<Point<dim>> &positions,
                                 const std::vector<unsigned int> &component_indices,
                                 std::vector<bool> &should_be_constrained,
                                 std::vector<double> &solution) override;

      private:
        std::shared_ptr<const ::WorldBuilder::World> world_builder;
    };
  }
}

#endif
