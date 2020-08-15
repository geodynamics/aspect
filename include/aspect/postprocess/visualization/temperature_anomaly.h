/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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
                                std::vector<Vector<double> > &computed_quantities) const override;

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
           */
          unsigned int n_slices;
          /**
           * Vector of temperature depth average values, padded to include two ghost values
           * above and below the surface and bottom of the domain to allow interpolation at
           * all depths.
           */
          std::vector<double> padded_temperature_depth_average;
          /**
           * Whether to extrapolate temperatures above/below the first/last depth-average slice
           * or, alternatively, interpolate above the center of the first slice using the surface
           * temperature or below the last slice using the bottom temperature.
           */
          bool extrapolate_surface;
          bool extrapolate_bottom;
      };
    }
  }
}

#endif
