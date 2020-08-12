/*
  Copyright (C) 2015 - 2019 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_visualization_stress_h
#define _aspect_postprocess_visualization_stress_h

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
       * A class derived from DataPostprocessor that takes an output vector
       * and computes a variable that represents the 3 or 6 independent
       * components (in 2d and 3d, respectively) of the stress tensor at every
       * point. The stress is defined as $2 \eta (\varepsilon(\mathbf u)
       * - \tfrac 13 \textrm{trace}\ \varepsilon(\mathbf u) \mathbf 1) +pI =
       * 2\eta (\varepsilon(\mathbf u) - \frac 13 (\nabla \cdot \mathbf u)
       * \mathbf I) + pI$.  The second term in the parentheses is zero if the
       * model is incompressible.
       *
       * The member functions are all implementations of those declared in the
       * base class. See there for their meaning.
       */
      template <int dim>
      class Stress
        : public DataPostprocessorTensor<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          /**
           * Constructor.
           */
          Stress ();

          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double> > &computed_quantities) const override;
      };
    }
  }
}

#endif
