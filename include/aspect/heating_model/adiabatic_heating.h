/*
  Copyright (C) 2014 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef __aspect__heating_model_adiabatic_heating_h
#define __aspect__heating_model_adiabatic_heating_h

#include <aspect/heating_model/interface.h>

namespace aspect
{
  namespace HeatingModel
  {
    using namespace dealii;

    /**
     * A class that implements a standard adiabatic heating rate.
     * This adds the term from adiabatic compression heating
     *    + alpha T (u . nabla p)
     * where we use the definition of
     *    alpha = - 1/rho drho/dT
     * Note: this term is often simplified using the relationship
     *    rho g = -nabla p
     * to yield
     *    - alpha rho T (u . g)
     * However, we do not use this simplification here, see the
     * comment in the manual in the section on the governing
     * equations.
     *
     * @ingroup HeatingModels
     */
    template <int dim>
    class AdiabaticHeating : public Interface<dim>
    {
      public:
        /**
         * Compute the heating model outputs for this class.
         */
        virtual
        void
        evaluate (const typename aspect::MaterialModel::Interface<dim>::MaterialModelInputs &material_model_inputs,
                  const typename aspect::MaterialModel::Interface<dim>::MaterialModelOutputs &material_model_outputs,
                  HeatingModel::HeatingModelOutputs &heating_model_outputs) const;

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
        virtual
        void
        parse_parameters (ParameterHandler &prm);
        /**
         * @}
         */

      private:
        bool simplified_adiabatic_heating;
    };
  }
}


#endif
