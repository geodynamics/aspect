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


#ifndef __aspect__heating_model_latent_heat_h
#define __aspect__heating_model_latent_heat_h

#include <aspect/heating_model/interface.h>

namespace aspect
{
  namespace HeatingModel
  {
    using namespace dealii;

    /**
     * A class that implements a standard formulation of latent heat.
     * This includes a left hand side and a right hand side term:
     *
     * The left hand side is
     *   $ -\rho T frac{\partial S}{\partial T} \frac{D T}{D t}$
     * and is added to the rho cp term.
     *
     * The right-hand side term from latent heating is
     *   $\frac{\partial S}{\partial p} T \rho (v \dot \nabla p)$.
     *
     * Formulation modified after Christensen & Yuen, 1985.
     *
     * @ingroup HeatingModels
     */
    template <int dim>
    class LatentHeat : public Interface<dim>
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
