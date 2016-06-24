/*
  Copyright (C) 2011-2015 by the authors of the ASPECT code.

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


#ifndef __aspect__postprocess_boundary_pressures_h
#define __aspect__postprocess_boundary_pressures_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes the laterally averaged pressure
     * at the top and bottom boundaries of the solution.  This is useful
     * for calculating the dyanamic topography at those surfaces.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class BoundaryPressures : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for the laterally averaged pressure at
         * the top and bottom of the domain.
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);

        /**
         * Query the pressure at the top surface
         */
        double pressure_at_top() const;

        /**
         * Query the pressure at the bottom surface
         */
        double pressure_at_bottom() const;

      private:

        /**
         * Pressure at the top of the domain.
         * Filled when execute() is called.
         */
        double top_pressure;

        /**
         * Pressure at the bottom of the domain.
         * Filled when execute() is called.
         */
        double bottom_pressure;
    };
  }
}


#endif
