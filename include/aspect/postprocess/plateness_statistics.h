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


#ifndef _aspect_postprocess_plateness_statistics_h
#define _aspect_postprocess_plateness_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {
    /**
     * A postprocessor that computes surface plateness diagnostics on the top
     * boundary using the second invariant of the deviatoric strain-rate
     * tensor. It computes F80 and F90, the fractions of the top-boundary area
     * required to account for 80% and 90% of the total surface deformation,
     * and the corresponding plateness values p80 and p90.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class PlatenessStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Compute the surface plateness statistics and add them to the
         * statistics object.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

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
         * The reference surface-area fraction used to compute plateness as
         * p = 1 - F/reference_fraction.
         */
        double reference_fraction;
    };
  }
}


#endif
