// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_postprocess_visualization_grain_lag_angle_h
#define _aspect_postprocess_visualization_grain_lag_angle_h

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {

      /**
       * This postprocessor calculates and outputs the angle between the ~infinite
       * strain axis and the velocity. Kaminski & Ribe (2002, Gcubed) call this quantity
       * $\Theta$ and define it as
       * $\Theta = \cos^{-1}(\hat{u}\cdot\hat{e})$
       * where $\hat{u}=\vec{u}/|{u}|$, $\vec{u}$ is the local flow velocity, and
       * $\hat{e}$ is the local infinite strain axis, which we calculate as the
       * first eigenvector of the "left stretch" tensor.
       * $\Theta$ can be used to calculate the grain orientation lag parameter (GOL).
       * Calculating GOL also requires the ISA rotation timescale ($\tau_{ISA}$).
       * GOL is not calculated within ASPECT
       * right now because it is proportional to the spatial gradient of $\Theta$, but in the
       * future that calculation could be implemented in a material model with CopyOutputs
       * (once they exist). By tracking $\Theta$ as a CopyOutput (ie, a compositional field
       * holding a calculated value that gets copied over instead of solved for), the spatial
       * gradient of $\Theta$ could be calculated for the previous timestep by obtaining the
       * old solution for the input material model. That gradient, and also the time derivative
       * of $\Theta$, could then be used to calculate GOL at the previous timestep; $\Theta$ could
       * be updated at the current timestep; and both quantities could be stored in CopyOutputs
       * to step forward in time. Basically, the calculation of GOL would have to lag one
       * timestep behind the other quantities in order to get the gradients, but we're
       * often interested in GOL in a steady-state flow anyway.
       */
      template <int dim>
      class GrainLagAngle: public CellDataVectorCreator<dim>, public SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          GrainLagAngle();

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
