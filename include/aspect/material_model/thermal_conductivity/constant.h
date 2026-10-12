// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2025-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_material_model_thermal_conductivity_constant_h
#define _aspect_material_model_thermal_conductivity_constant_h

#include <aspect/material_model/thermal_conductivity/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {
      /**
       * A class that implements a constant thermal conductivity.
       *
       * @ingroup MaterialModels
       */
      template <int dim>
      class Constant : public Interface<dim>
      {
        public:
          /**
           * Function to compute the thermal conductivities in @p out given the
           * inputs in @p in.
           */
          void evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                         MaterialModel::MaterialModelOutputs<dim> &out) const override;

          /**
           * Declare the parameters this plugin takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm) override;

        private:
          /**
           * The thermal conductivity.
           *
           * This variable is read from the parameter file through a parameter called 'Thermal conductivity'.
           */
          double k;
      };
    }
  }
}

#endif
