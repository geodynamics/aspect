/*
  Copyright (C) 2022 by the authors of the ASPECT code.

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

// create a visualization postprocessor that requires something else
//
// this works just like the plugin_dependency* tests, just that it's a
// viz postprocessor, not a regular postprocessor that has
// dependencies

#include "../benchmarks/solcx/solcx.cc"
#include <aspect/postprocess/visualization.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      class MyPostprocessor
        : public DataPostprocessorScalar<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          MyPostprocessor ()
            :
            DataPostprocessorScalar<dim> ("my_postprocessor",
                                          update_default)
          {}

          virtual
          void
          evaluate_vector_field (const DataPostprocessorInputs::Vector<dim> &input_data,
                                 std::vector<Vector<double>>               &computed_quantities) const
          {
            Assert (computed_quantities[0].size() == 1, ExcInternalError());

            for (unsigned int q=0; q<computed_quantities.size(); ++q)
              {
                // not important what we do here :-)
                computed_quantities[q](0) = 0;
              }
          }

          virtual
          std::list<std::string>
          required_other_postprocessors () const
          {
            // select a postprocessor that is not selected in the .prm file
            std::list<std::string> deps;
            deps.push_back ("velocity statistics");
            return deps;
          }
      };
    }
  }
}


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(MyPostprocessor,
                                                  "my postprocessor",
                                                  ".")
    }
  }
}
