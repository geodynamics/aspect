/*
  Copyright (C) 2017 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_mass_error_statistics_h
#define _aspect_postprocess_mass_error_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {
    /**
     * A postprocessor that computes some statistics about the "
     * error of the mass conservation equation. More precisely, it
     * computes the norm and the integral of the equation
     * $\\rho \\Nabla \\mathbf u + \\Nabla \\rho \\cdot \\mathbf u$.
     * Therefore, its results indicate how accurate the chosen approximation
     * for compressibility is for the given model setup.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class MassErrorStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for some velocity statistics.
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);
    };
  }
}


#endif
