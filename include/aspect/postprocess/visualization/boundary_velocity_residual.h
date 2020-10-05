/*
  Copyright (C) 2020 by the authors of the ASPECT code.

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
      using namespace dealii;
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
           * @copydoc DataPostprocessorVector<dim>::evaluate_vector_field()
           */
          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double> > &computed_quantities) const override;

          /**
           * Let the postprocessor manager know about the other postprocessors
           * this one depends on. Specifically, the boundary velocity residual statistics postprocessor.
           */
          std::list<std::string>
          required_other_postprocessors() const override;
      };
    }
  }
}

#endif