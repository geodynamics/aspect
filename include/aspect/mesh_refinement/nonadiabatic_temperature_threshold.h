// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2015-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


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
         *
         * This variable is read from the parameter file through a parameter called 'Threshold'.
         */
        double threshold;

        /**
         * What type of temperature anomaly should be considered when
         * evaluating against the threshold: Only negative anomalies
         * (subadiabatic temperatures), only positive anomalies
         * (superadiabatic temperatures) or the absolute value of the
         * nonadiabatic temperature.
         *
         * This variable is read from the parameter file through a parameter called 'Temperature anomaly type'.
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
