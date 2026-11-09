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


#ifndef _aspect_postprocess_visualization_boundary_strain_rate_residual_h
#define _aspect_postprocess_visualization_boundary_strain_rate_residual_h

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
       * and computes a variable that represents the residual between
       * the a reference strain rate dataset and the strain rate computed for the model domain.
       * We only compute the residual between the second invariant of the strain rate
       * and therefore, the input data is an ascii data file with columns corresponding
       * to coordinates and one data column representing the second invariant of the strain rate.
       *
       * This quantity only makes sense at the surface of the domain.
       * Thus, the value is set to zero in all the cells inside of the domain.
       *
       * The member functions are all implementations of those declared in the
       * base class. See there for their meaning.
       */
      template <int dim>
      class BoundaryStrainRateResidual
        : public DataPostprocessorScalar<dim>,
          public SurfaceOnlyVisualization<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          BoundaryStrainRateResidual ();

          /**
           * Evaluate the strain rate residual for each component at the top face of the current cell
           * if it lies at the surface boundary.
           *
           * @copydoc dealii::DataPostprocessor<dim>::evaluate_vector_field()
           */
          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double>> &computed_quantities) const override;

          /**
           * Let the postprocessor manager know about the other postprocessors
           * this one depends on. Specifically, the boundary strain rate residual statistics postprocessor.
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
