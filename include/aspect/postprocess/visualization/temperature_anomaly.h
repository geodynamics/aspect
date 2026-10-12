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


#ifndef _aspect_postprocess_visualization_temperature_anomaly_h
#define _aspect_postprocess_visualization_temperature_anomaly_h

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>

#include <deal.II/numerics/data_postprocessor.h>


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
      class TemperatureAnomaly
        : public DataPostprocessorScalar<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          TemperatureAnomaly ();

          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double>> &computed_quantities) const override;

          void
          update () override;

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
           * Number of slices to use when computing depth average of temperature.
           *
           * This variable is read from the parameter file through a parameter called 'Number of depth slices'.
           */
          unsigned int n_slices;

          /**
           * Whether to extrapolate temperatures above/below the first/last depth-average slice
           * or, alternatively, interpolate above the center of the first slice using the surface
           * temperature or below the last slice using the bottom temperature.
           * This variable is read from the parameter file through a parameter called 'Use minimal temperature for surface'.
           */
          bool extrapolate_surface;

          /**
           * This variable is read from the parameter file through a parameter called 'Use maximal temperature for bottom'.
           */
          bool extrapolate_bottom;

          /**
           * Vector of temperature depth average values, padded to include two ghost values
           * above and below the surface and bottom of the domain to allow interpolation at
           * all depths.
           *
           * This object is computed in update() and then used on every cell
           * where evaluate_vector_field() is called.
           */
          std::vector<double> padded_temperature_depth_average;
      };
    }
  }
}

#endif
