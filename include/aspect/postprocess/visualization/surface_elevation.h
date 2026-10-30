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


#ifndef _aspect_postprocess_visualization_surface_elevation_h
#define _aspect_postprocess_visualization_surface_elevation_h

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * A class derived from DataPostprocessorScalar that computes a
       * variable that represents the elevation of points on the ``top''
       * surface relative to the geometry's undeformed state, and outputs
       * it on boundary faces. On boundary faces not marked as ``top'',
       * this class outputs a zero elevation.
       *
       * The member functions are all implementations of those declared in the
       * base class. See there for their meaning.
       */
      template <int dim>
      class SurfaceElevation
        : public DataPostprocessorScalar<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>,
          public SurfaceOnlyVisualization<dim>
      {
        public:
          SurfaceElevation();

          /**
           * Evaluate the dynamic topography for the current cell.
           *
           * @copydoc dealii::DataPostprocessor<dim>::evaluate_vector_field()
           */
          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double>> &computed_quantities) const override;
      };
    }
  }
}

#endif
