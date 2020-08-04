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


#ifndef _aspect_postprocess_visualization_surface_stress_h
#define _aspect_postprocess_visualization_surface_stress_h

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
       * point of the surface of the domain.
       * The stress is defined as $2 \eta (\varepsilon(\mathbf u)
       * - \tfrac 13 \textrm{trace}\ \varepsilon(\mathbf u) \mathbf 1) +pI =
       * 2\eta (\varepsilon(\mathbf u) - \frac 13 (\nabla \cdot \mathbf u)
       * \mathbf I) + pI$.  The second term in the parentheses is zero if the
       * model is incompressible.
       *
       * The member functions are all implementations of those declared in the
       * base class. See there for their meaning.
       */
      template <int dim>
      class SurfaceStress
        : public DataPostprocessor<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>,
          public SurfaceOnlyVisualization<dim>
      {
        public:
          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double> > &computed_quantities) const override;

          /**
           * Return the vector of strings describing the names of the computed
           * quantities. Given the purpose of this class, this is a vector
           * with entries all equal to the name of the plugin.
           */
          std::vector<std::string> get_names () const override;

          /**
           * This functions returns information about how the individual
           * components of output files that consist of more than one data set
           * are to be interpreted. The returned value is
           * DataComponentInterpretation::component_is_scalar repeated
           * SymmetricTensor::n_independent_components times. (These
           * components should really be part of a symmetric tensor, but
           * deal.II does not allow marking components as such.)
           */
          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          get_data_component_interpretation () const override;

          /**
           * Return which data has to be provided to compute the derived
           * quantities. The flags returned here are the ones passed to the
           * constructor of this class.
           */
          UpdateFlags get_needed_update_flags () const override;
      };
    }
  }
}

#endif
