/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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
           */
          unsigned int n_slices;


      };
    }
  }
}

#endif
