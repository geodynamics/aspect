// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_postprocess_visualization_density_anomaly_h
#define _aspect_postprocess_visualization_density_anomaly_h

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * A class derived from DataPostprocessor that outputs the temperature
       * anomaly (temperature-depth average of temperature).
       */
      template <int dim>
      class DensityAnomaly
        :
        public CellDataVectorCreator<dim>,
        public SimulatorAccess<dim>
      {
        public:
          DensityAnomaly ();

          /**
           * @copydoc CellDataVectorCreator<dim>::execute()
           */
          std::pair<std::string, std::unique_ptr<Vector<float>>>
          execute () const override;

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
           * Scheme chosen to define the average density as a
           * function of depth. Reference profile evaluates the material model
           * using the P-T profile defined by the reference adiabatic
           * conditions and the lateral average option calculates the average
           * density within a number n_slices of depth slices.
           */
          enum DensityScheme
          {
            reference_profile,
            lateral_average
          } average_density_scheme;

          /**
           * Number of slices to use when computing depth average of temperature.
           *
           * This variable is read from the parameter file through a parameter called 'Number of depth slices'.
           */
          unsigned int n_slices;


      };
    }
  }
}

#endif
