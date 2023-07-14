/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_composition_velocity_statistics_h
#define _aspect_postprocess_composition_velocity_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes the root mean square velocity
     * over the area spanned by each of the compositional fields.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class CompositionVelocityStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution to compute the RMS velocity for each composition.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */
        /**
         * Declare the parameters this class takes through input files.
         */
        static void
        declare_parameters(ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters(ParameterHandler &prm) override;
        /**
         * @}
         */

      private:
        std::vector<std::string> selected_fields;
    };
  }
}


#endif
