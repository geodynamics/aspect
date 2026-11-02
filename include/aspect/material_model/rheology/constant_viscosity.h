// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2019-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_material_model_rheology_constant_viscosity_h
#define _aspect_material_model_rheology_constant_viscosity_h

#include <aspect/global.h>
#include <deal.II/base/parameter_handler.h>

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
           *
           * This variable is read from the parameter file through a parameter called 'Viscosity'.
           */
          double viscosity;
      };
    }
  }
}
#endif
