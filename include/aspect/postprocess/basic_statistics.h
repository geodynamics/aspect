/*
  Copyright (C) 2011 - 2017 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_basic_statistics_h
#define _aspect_postprocess_basic_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that outputs some simplified statistics
     * like the Rayleigh number and other quantities that only
     * make sense in certain model setups. The output is written
     * after completing initial adaptive refinement steps.
     * The postprocessor assumes a point at the surface at the adiabatic
     * surface temperature and pressure is a reasonable reference condition
     * for computing these properties. Furthermore, the Rayleigh
     * number is computed using the model depth (i.e. not the
     * radius of the Earth), as we need a definition that is geometry
     * independent. Take care when comparing these values to published
     * studies and make sure they use the same definitions.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class BasicStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for some general statistics.
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);
    };
  }
}


#endif
