/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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


/**
 * This postprocessor calculates and outputs two quantities:

 *   -> theta, the angle between the ~infinite strain axis and the velocity
 *   -> tauISA, the timescale for the rotation of grains toward the infinite strain axis.

 * Those two quantites can be used to calculate the grain orientation lag parameter (GOL)
 * defined by Kaminski and Ribe (2002, Gcubed). GOL is not calculated within ASPECT
 * right now because it is proportional to the spatial gradient of theta, but in the
 * future that calculation could be implemented in a material model with CopyOutputs
 * (once they exist). By tracking theta as a CopyOutput (ie, a compositional field
 * holding a calculated value that gets copied over instead of solved for), the spatial
 * gradient of theta could be calculated for the previous timestep by obtaining the
 * old solution for the input material model. That gradient, and also the time derivative
 * of theta, could then be used to calculate GOL at the previous timestep; theta could
 * be updated at the current timestep; and both quantities could be stored in CopyOutputs
 * to step forward in time. Basically, the calculation of GOL would have to lag one
 * timestep behind the other quantities in order to get the gradients, but we're
 * often interested in GOL in a steady-state flow anyway.
*/


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * A derived class that implements a function that provides the computed
       * timescale for the rotation of grains toward the infinite strain axis
       * for graphical output.
       */
      template<int dim>
      class ISARotationTimescale: public CellDataVectorCreator<dim>, public SimulatorAccess<dim>
      {
        public:

          virtual std::pair<std::string, Vector<float> *>
          execute() const;

      };
    }
  }
}

#endif
