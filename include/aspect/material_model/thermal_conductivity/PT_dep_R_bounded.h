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

#ifndef _aspect_material_model_thermal_conductivity_PT_dep_R_bounded_h
#define _aspect_material_model_thermal_conductivity_PT_dep_R_bounded_h

#include <aspect/material_model/thermal_conductivity/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {
      using namespace dealii;

      /**
       * A base class for parameterizations of material models. Classes derived
       * from this class will need to implement functions that provide material
       * parameters such as the viscosity, density, etc, typically as a function
       * of position, temperature and pressure at that location.
       *
       * Implementing a material model requires you to override evaluate() and fill the output
       * argument struct instead of implementing the functions viscosity(),
       * density(), etc.. In this case, all other functions are being ignored.
       *
       * In all cases, model_dependence values, is_compressible()
       * need to be implemented.
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
           */
          double k;
      };
    }
  }
}

#endif
