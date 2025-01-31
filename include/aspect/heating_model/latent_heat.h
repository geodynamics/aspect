/*
  Copyright (C) 2011 - 2021 by the authors of the ASPECT code.

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


#ifndef _aspect_heating_model_latent_heat_h
#define _aspect_heating_model_latent_heat_h

#include <aspect/heating_model/interface.h>

namespace aspect
{
  namespace HeatingModel
  {
    /**
     * A class that implements a standard formulation of latent heat.
     * This includes a left hand side and a right hand side term:
     *
     * The left hand side is
     *   $ -\rho T \frac{\partial S}{\partial T} \frac{D T}{D t}$
     * so that we can add
     *   $ -\rho T \frac{\partial S}{\partial T} $
     * to the $\rho C_p$ term.
     *
     * The right-hand side term from latent heating is
     *   $\frac{\partial S}{\partial p} T \rho (u \dot \nabla p)$.
     *
     * T, u, and p are the solutions from the previous time step or
     * are extrapolated from there, depending on what is provided
     * in the input arguments of this function.
     *
     * Formulation modified after Christensen, Ulrich R. & Yuen,
     * David A.: Layered convection induced by phase transitions,
     * Journal of Geophysical Research: Solid Earth (1985).
     *
     * Also see the Equations section in the manual.
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
        void
        evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                  const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                  HeatingModel::HeatingModelOutputs &heating_model_outputs) const override;
    };
  }
}


#endif
