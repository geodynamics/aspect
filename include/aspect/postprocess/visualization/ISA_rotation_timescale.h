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


#ifndef _aspect_postprocess_visualization_ISA_rotation_timescale_h
#define _aspect_postprocess_visualization_ISA_rotation_timescale_h

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>





namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * This postprocessor calculates and outputs the timescale for the rotation
       * of grains toward the infinite strain axis. Kaminski and Ribe (2002, Gcubed)
       * call this quantity $\tau_{ISA}$, and define it as
       * $\tau_{ISA} \approx \frac{1}{\dot{\epsilon}}$
       * where $\dot{\epsilon}$ is the largest eigenvalue of the strain rate tensor.
       * It can be used, along with the grain lag angle
       * ($\Theta$), to calculate the grain orientation lag parameter (GOL).
       * GOL is not calculated within ASPECT
       * right now because it is proportional to the spatial gradient of theta, but in the
       * future that calculation could be implemented in a material model with CopyOutputs
       * (once they exist). For more thoughts on that, see the documentation for the
       * grain lag angle postprocessor.
       */
      template<int dim>
      class ISARotationTimescale: public CellDataVectorCreator<dim>, public SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          ISARotationTimescale();

          /**
           * @copydoc CellDataVectorCreator<dim>::execute()
           */
          std::pair<std::string, Vector<float> *>
          execute() const override;

      };
    }
  }
}

#endif
