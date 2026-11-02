// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_postprocess_visualization_depth_including_mesh_deformation_h
#define _aspect_postprocess_visualization_depth_including_mesh_deformation_h

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
       * A class derived from DataPostprocessor that outputs the depth
       * including mesh deformation according to the geometry model.
       */
      template <int dim>
      class DepthIncludingMeshDeformation
        : public DataPostprocessorScalar<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          DepthIncludingMeshDeformation ();

          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double>> &computed_quantities) const override;

          /**
           * Let the postprocessor manager know about the other postprocessors
           * this one depends on. Specifically, the Geoid postprocessor.
           */
          std::list<std::string>
          required_other_postprocessors() const override;
      };
    }
  }
}

#endif
