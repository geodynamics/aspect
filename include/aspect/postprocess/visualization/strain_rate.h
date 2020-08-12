/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_visualization_strain_rate_h
#define _aspect_postprocess_visualization_strain_rate_h

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
       * and computes a variable that represents the strain rate at every
       * point. The scalar strain rate is defined as $\sqrt{ (\varepsilon -
       * \tfrac 13 \textrm{trace}\ \varepsilon \mathbf 1) : \varepsilon -
       * \tfrac 13 \textrm{trace}\ \varepsilon \mathbf 1}$.
       *
       * The member functions are all implementations of those declared in the
       * base class. See there for their meaning.
       */
      template <int dim>
      class StrainRate
        : public DataPostprocessorScalar<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          /**
           * Constructor.
           */
          StrainRate ();

          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double> > &computed_quantities) const override;
      };
    }
  }
}

#endif
