/*
  Copyright (C) 2017 - 2021 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_visualization_spd_factor_h
#define _aspect_postprocess_visualization_spd_factor_h

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
       * and computes a variable that represents the spd factor at every point.
       * The spd factor is a factor which scales a part of the Jacobian used
       * for the Newton solver to make sure that the Jacobian remains positive
       * definite.
       *
       * The member functions are all implementations of those declared in the
       * base class. See there for their meaning.
       */
      template <int dim>
      class SPD_Factor
        : public DataPostprocessorScalar<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          SPD_Factor ();

          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double>> &computed_quantities) const override;

          /**
           * Read the parameters this class declares from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm) override;

      };
    }
  }
}

#endif
