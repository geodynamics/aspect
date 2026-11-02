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

#ifndef _aspect_material_model_simpler_h
#define _aspect_material_model_simpler_h

#include <aspect/material_model/interface.h>
#include <aspect/material_model/rheology/constant_viscosity.h>
#include <aspect/material_model/equation_of_state/linearized_incompressible.h>
#include <aspect/material_model/thermal_conductivity/constant.h>

namespace aspect
{
  namespace MaterialModel
  {
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
        ThermalConductivity::Constant<dim> thermal_conductivity;
        Rheology::ConstantViscosity constant_rheology;
        EquationOfState::LinearizedIncompressible<dim> equation_of_state;
    };

  }
}

#endif
