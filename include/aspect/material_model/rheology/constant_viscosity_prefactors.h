// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2020-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_material_model_rheology_constant_viscosity_prefactors_h
#define _aspect_material_model_rheology_constant_viscosity_prefactors_h

#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      /**
       * A class that handles multiplication of viscosity for a given compositional
       * field. The multiplication factors for each composition (constant viscosity
       * prefactors) are also declared and parsed in this class.
       */
      template <int dim>
      class ConstantViscosityPrefactors : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          ConstantViscosityPrefactors();

          /**
           * Declare the parameters this function takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm);

          /**
           * Compute the viscosity.
           */
          double
          compute_viscosity (const double base_viscosity,
                             const unsigned int composition_index) const;

        private:
          /**
           * The constant viscosity prefactors, which are read in
           * from the input file by the parse_parameters() function.
           * The total number of prefactors will be equal to one
           * plus the number of compositional fields. The prefactor
           * for a given compositional field is multiplied with a
           * base_viscosity value provided by the material model, which
           * is then returned to the material model.
           *
           * This variable is read from the parameter file through a parameter called 'Constant viscosity prefactors'.
           */
          std::vector<double> constant_viscosity_prefactors;
      };
    }
  }
}
#endif
