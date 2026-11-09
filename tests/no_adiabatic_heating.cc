// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2022-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

// use the same postprocessing facilities as for the
// 'compressibility_iterated_stokes' testcase
#include "compressibility_iterated_stokes.cc"

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    class NoAdiabaticHeating : public CompressibilityIteratedStokes<dim>
    {
      public:
        virtual double thermal_expansion_coefficient (const double,
                                                      const double,
                                                      const std::vector<double> &,
                                                      const Point<dim> &) const
        {
          return 0;
        }
    };

  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(NoAdiabaticHeating,
                                   "no adiabatic heating",
                                   "As described in the .prm file.")
  }
}
