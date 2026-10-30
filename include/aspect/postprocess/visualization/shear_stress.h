// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2015-2021 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_postprocess_visualization_shear_stress_h
#define _aspect_postprocess_visualization_shear_stress_h

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
       * components (in 2d and 3d, respectively) of the shear stress tensor at
       * every point. The shear stress is defined as $2 \eta
       * (\varepsilon(\mathbf u) - \tfrac 13 \textrm{trace}\
       * \varepsilon(\mathbf u) \mathbf 1) = 2\eta (\varepsilon(\mathbf u) -
       * \frac 13 (\nabla \cdot \mathbf u) \mathbf I)$ and differs from the
       * full stress by the absence of the pressure.  The second term in the
       * difference is zero if the model is incompressible.
       * If elasticity is used, the elastic contribution is being accounted for.
       *
       * The member functions are all implementations of those declared in the
       * base class. See there for their meaning.
       */
      template <int dim>
      class ShearStress
        : public DataPostprocessorTensor<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          /**
           * Constructor.
           */
          ShearStress ();

          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double>> &computed_quantities) const override;
      };
    }
  }
}

#endif
