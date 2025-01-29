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
