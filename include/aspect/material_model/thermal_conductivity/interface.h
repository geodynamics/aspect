// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_material_model_thermal_conductivity_interface_h
#define _aspect_material_model_thermal_conductivity_interface_h

#include <aspect/global.h>
#include <aspect/plugins.h>
#include <aspect/material_model/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {
      /**
       * A base class for parametrizations of the thermal conductivity. Classes derived
       * from this class will need to implement a function that computes the thermal
       * conductivities in @p out given the inputs in @p in. Derived classes can in addition
       * implement the functions of the base class Plugins::InterfaceBase as needed (e.g.
       * to read in input parameters or update the parametrization at the beginning of time steps).
       *
       * @ingroup MaterialModels
       */
      template <int dim>
      class Interface : public Plugins::InterfaceBase
      {
        public:
          /**
           * Function to compute the thermal conductivities in @p out given the
           * inputs in @p in.
           */
          virtual
          void evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                         MaterialModel::MaterialModelOutputs<dim> &out) const = 0;
      };
    }
  }
}

#endif
