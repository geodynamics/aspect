/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_simpler_h
#define _aspect_material_model_simpler_h

#include <aspect/material_model/interface.h>
#include <aspect/material_model/rheology/constant_viscosity.h>
#include <aspect/material_model/equation_of_state/linearized_incompressible.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model that consists of globally constant values for all
     * material parameters except the density, which depends linearly on the
     * temperature. The model is considered incompressible.
     *
     * This material model implements what the "Simple" model was originally
     * intended to do, before it got too complicated.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class Simpler : public Interface<dim>
    {
      public:

        bool is_compressible () const override;

        double reference_viscosity () const override;

        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override;


        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */
        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * @}
         */

      private:
        double k_value;

        Rheology::ConstantViscosity constant_rheology;
        EquationOfState::LinearizedIncompressible<dim> equation_of_state;
    };

  }
}

#endif
