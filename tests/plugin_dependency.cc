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

// create a postprocessor that does exactly the same as
// SolCxPostprocessor, but also requires something else

#include "../benchmarks/solcx/solcx.cc"

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    class NewPostprocessor : public aspect::InclusionBenchmark::SolCxPostprocessor<dim>
    {
      public:
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


namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(NewPostprocessor,
                                  "new postprocessor",
                                  ".")
  }
}
