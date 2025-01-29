/*
  Copyright (C) 2019 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_rheology_constant_viscosity_h
#define _aspect_material_model_rheology_constant_viscosity_h

#include <aspect/global.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      class ConstantViscosity
      {
        public:
          /**
           * Constructor. Initializes viscosity to NaN.
           */
          ConstantViscosity();

          /**
           * Declare the parameters this function takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm,
                              const double default_viscosity = 1e21);

          /**
           * Read the parameters from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm);

          /**
           * Compute the viscosity, in this case just a constant.
           */
          double
          compute_viscosity () const;

        private:
          /**
           * The constant viscosity that defines this rheology. It
           * is read from the input file by the parse_parameters()
           * function.
           */
          double viscosity;
      };
    }
  }
}
#endif
