// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


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

        /**
         * Specify which material model outputs the heating model requires
         * for computing the heating terms.
         */
        MaterialModel::MaterialProperties::Property
        get_required_properties () const override;
    };
  }
}


#endif
