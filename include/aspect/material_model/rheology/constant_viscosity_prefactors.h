/*
  Copyright (C) 2020 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_rheology_constant_viscosity_prefactors_h
#define _aspect_material_model_rheology_constant_viscosity_prefactors_h

#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

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
           */
          std::vector<double> constant_viscosity_prefactors;
      };
    }
  }
}
#endif
