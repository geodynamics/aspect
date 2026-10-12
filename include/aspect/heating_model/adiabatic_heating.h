// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_heating_model_adiabatic_heating_h
#define _aspect_heating_model_adiabatic_heating_h

#include <aspect/heating_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace HeatingModel
  {
    /**
     * A class that implements a standard adiabatic heating rate.
     *
     * This adds the term from adiabatic compression heating
     *    $ \alpha T (\mathbf u \cdot \nabla p) $
     * where we use the definition of
     *    $ \alpha = - \frac{1}{\rho} \frac{\partial \rho}{\partial T} $
     * Note: this term is often simplified using the relationship
     *    $ \rho \mathbf g = - \nabla p $
     * to yield
     *    $ - \alpha \rho T (\mathbf u \cdot \mathbf g) $
     *
     * The user can specify if the simplification should be used
     * by setting the corresponding flag in the input file.
     * Also see the Equations section in the manual.
     *
     * @ingroup HeatingModels
     */
    template <int dim>
    class AdiabaticHeating : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Returns whether we compute a simplified adiabatic heating term by
         * neglecting dynamic pressure effects. I.e. we use
         *     $ \alpha T (\mathbf u \cdot \nabla p) $
         * as adiabatic heating term if this function returns false, and
         *    $ - \alpha \rho T (\mathbf u \cdot \mathbf g) $
         * if it returns true.
         */
        bool
        use_simplified_adiabatic_heating() const;

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
        /**
         * This variable is read from the parameter file through a parameter called 'Use simplified adiabatic heating'.
         */
        bool simplified_adiabatic_heating;
    };
  }
}


#endif
