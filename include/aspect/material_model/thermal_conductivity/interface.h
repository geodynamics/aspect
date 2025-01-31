/*
  Copyright (C) 2025 - by the authors of the ASPECT code.

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
