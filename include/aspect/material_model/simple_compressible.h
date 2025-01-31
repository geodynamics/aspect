/*
  Copyright (C) 2014 - 2022 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_simple_compressible_h
#define _aspect_material_model_simple_compressible_h

#include <aspect/material_model/interface.h>
#include <aspect/material_model/rheology/constant_viscosity.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * A material model that consists of globally constant values for the
     * viscosity, thermal conductivity, thermal expansivity
     * and compressibility. The density decays linearly with the
     * temperature and increases exponentially with pressure.
     *
     * The formulation for the density assumes that the compressibility
     * provided by the user is the adiabatic compressibility ($\beta_S$).
     * The thermal expansivity and isentropic compressibility implied by
     * the pressure and temperature dependence are equal to the
     * user-defined constant values only along the reference isentrope, and
     * there is also an implicit pressure dependence to the heat capacity
     * $C_p$ via Maxwell's relations.
     *
     * The model is considered incompressible or compressible, depending on
     * the compressibility.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class SimpleCompressible : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate material properties.
         */
        void evaluate(const MaterialModelInputs<dim> &in,
                      MaterialModelOutputs<dim> &out) const override;

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the continuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        bool is_compressible () const override;
        /**
         * @}
         */

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
         * The reference density
         */
        double reference_rho;

        /**
         * The constant thermal expansivity
         */
        double thermal_alpha;

        /**
         * The constant specific heat
         */
        double reference_specific_heat;

        /**
         * The constant compressibility.
         */
        double reference_compressibility;

        /**
         * The constant thermal conductivity.
         */
        double k_value;

        /**
         * The constant rheology model.
         */
        Rheology::ConstantViscosity constant_rheology;
    };

  }
}

#endif
