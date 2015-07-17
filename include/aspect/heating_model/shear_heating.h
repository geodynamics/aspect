/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
     * Add the term
     *    $  2 \eta \left( \varepsilon - \frac{1}{3} \text{tr}
     *       \varepsilon \mathbf 1 \right) : \left( \varepsilon - \frac{1}{3}
     *       \text{tr} \varepsilon \mathbf 1 \right) $
     *
     * Also see the Equations section in the manual.
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
        evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                  const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
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
