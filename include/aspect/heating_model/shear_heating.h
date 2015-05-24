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


#ifndef __aspect__heating_model_shear_heating_h
#define __aspect__heating_model_shear_heating_h

#include <aspect/heating_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace HeatingModel
  {
    using namespace dealii;

    /**
     * A class that implements a standard shear heating rate.
     *
     * Add the term 2*eta*(eps - 1/3*(tr eps)1):(eps - 1/3*(tr eps)1)
     *
     * we can multiply this out to obtain
     *      2*eta*(eps:eps - 1/3*(tr eps)^2)
     * and can then use that in the compressible case we have
     *      tr eps = div u
     *             = -1/rho u . grad rho
     * and by the usual approximation we make,
     *      tr eps = -1/rho drho/dp u . grad p
     *             = -1/rho drho/dp rho (u . g)
     *             = - drho/dp (u . g)
     *             = - compressibility rho (u . g)
     * to yield the final form of the term:
     *      2*eta [eps:eps - 1/3 (compressibility * rho * (u.g))^2]
     *
     * @ingroup HeatingModels
     */
    template <int dim>
    class ShearHeating : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
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

    };
  }
}


#endif
