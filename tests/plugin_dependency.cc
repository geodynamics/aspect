// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2022 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

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
