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


#ifndef _aspect_postprocess_visualization_seismic_anomalies_h
#define _aspect_postprocess_visualization_seismic_anomalies_h

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * A class derived that implements a function that provides the computed
       * computed seismic anomaly in $V_s$ for graphical output.
       */
      template <int dim>
      class SeismicVsAnomaly
        :
        public CellDataVectorCreator<dim>,
        public SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          SeismicVsAnomaly();

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
           * Scheme chosen to define the average seismic velocity as a
           * function of depth. Reference profile evaluates the material model
           * using the P-T profile defined by the reference adiabatic
           * conditions and the lateral average option calculates the average
           * velocity within a number n_slices of depth slices.
           */
          enum VelocityScheme
          {
            reference_profile,
            lateral_average
          } average_velocity_scheme;

          /**
           * Number of depth slices used to define average
           * seismic shear wave velocities from which anomalies
           * are calculated.
           *
           * This variable is read from the parameter file through a parameter called 'Number of depth slices'.
           */
          unsigned int n_slices;
      };



      /**
       * A class derived that implements a function that provides the computed
       * computed seismic anomaly in $V_p$ for graphical output.
       */
      template <int dim>
      class SeismicVpAnomaly
        :
        public CellDataVectorCreator<dim>,
        public SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          SeismicVpAnomaly();

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
           * Scheme chosen to define the average seismic velocity as a
           * function of depth. Reference profile evaluates the material model
           * using the P-T profile defined by the reference adiabatic
           * conditions and the lateral average option calculates the average
           * velocity within a number n_slices of depth slices.
           */
          enum VelocityScheme
          {
            reference_profile,
            lateral_average
          } average_velocity_scheme;

          /**
           * Number of depth slices used to define average
           * seismic compressional wave velocities from which anomalies
           * are calculated.
           *
           * This variable is read from the parameter file through a parameter called 'Number of depth slices'.
           */
          unsigned int n_slices;
      };
    }
  }
}

#endif
