// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2023 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_postprocess_visualization_artificial_viscosity_h
#define _aspect_postprocess_visualization_artificial_viscosity_h

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * A class derived that implements a function that provides the
       * artificial viscosity for the temperature equation on each cell for
       * graphical output.
       */
      template <int dim>
      class ArtificialViscosity : public CellDataVectorCreator<dim>,
        public SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          ArtificialViscosity();

          /**
           * @copydoc CellDataVectorCreator<dim>::execute()
           */
          std::pair<std::string, std::unique_ptr<Vector<float>>>
          execute () const override;
      };
    }
  }
}

#endif
