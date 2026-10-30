// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2024 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


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
      template <int dim>
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
          std::pair<std::string, std::unique_ptr<Vector<float>>>
          execute() const override;

      };
    }
  }
}

#endif
