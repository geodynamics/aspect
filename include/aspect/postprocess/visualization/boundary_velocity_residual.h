// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2020-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_postprocess_visualization_boundary_velocity_residual_h
#define _aspect_postprocess_visualization_boundary_velocity_residual_h

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

#include <deal.II/numerics/data_postprocessor.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * A class derived from DataPostprocessor that takes an output vector
       * and computes a variable that represents the velocity residual between
       * the input velocity data and the velocity computed for the model domain.
       * The input velocity could either be an ascii data file or a file generated
       * by  the GPlates model.
       * This quantity only makes sense at the surface of the domain.
       * Thus, the value is set to zero in all the cells inside of the domain.
       *
       * The member functions are all implementations of those declared in the
       * base class. See there for their meaning.
       */
      template <int dim>
      class BoundaryVelocityResidual
        : public DataPostprocessorVector<dim>,
          public SurfaceOnlyVisualization<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          BoundaryVelocityResidual ();

          /**
           * Evaluate the velocity residual for the current cell.
           *
           * @copydoc dealii::DataPostprocessor<dim>::evaluate_vector_field()
           */
          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double>> &computed_quantities) const override;

          /**
           * Let the postprocessor manager know about the other postprocessors
           * this one depends on. Specifically, the boundary velocity residual statistics postprocessor.
           */
          std::list<std::string>
          required_other_postprocessors() const override;

          std::string
          get_physical_units () const override;
      };
    }
  }
}

#endif
