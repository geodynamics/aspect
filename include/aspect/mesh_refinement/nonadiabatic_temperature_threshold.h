/*
  Copyright (C) 2015 - 2023 by the authors of the ASPECT code.

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


#ifndef _aspect_mesh_refinement_nonadiabatic_temperature_threshold_h
#define _aspect_mesh_refinement_nonadiabatic_temperature_threshold_h

#include <aspect/mesh_refinement/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MeshRefinement
  {

    /**
     * A class that implements a mesh refinement criterion based on the
     * nonadiabatic temperature. Every cell that contains a value
     * exceeding a threshold given in the input file is marked for
     * refinement.
     *
     * @ingroup MeshRefinement
     */
    template <int dim>
    class NonadiabaticTemperatureThreshold : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * After cells have been marked for coarsening/refinement, apply
         * additional criteria independent of the error estimate.
         */
        void
        tag_additional_cells () const override;

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
         * The thresholds that should be used for the nonadiabatic
         * temperature.
         */
        double threshold;

        /**
         * What type of temperature anomaly should be considered when
         * evaluating against the threshold: Only negative anomalies
         * (subadiabatic temperatures), only positive anomalies
         * (superadiabatic temperatures) or the absolute value of the
         * nonadiabatic temperature.
         */
        enum anomaly
        {
          negative_only,
          positive_only,
          absolute_value
        } temperature_anomaly_type;
    };
  }
}

#endif
