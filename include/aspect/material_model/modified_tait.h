/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_modified_tait_h
#define _aspect_material_model_modified_tait_h

#include <deal.II/base/function_lib.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * A compressible material model that implements the thermal modified Tait
     * equation of state as written in the paper of Holland and Powell, 2011
     * "An improved and extended internally consistent thermodynamic dataset
     * for phases of petrological interest, involving a new equation of state
     * for solids".
     *
     * Constant values are used for the thermal conductivity and viscosity.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class ModifiedTait : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
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
         * This model is compressible.
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
         * The Tait parameters
         */
        double tait_a, tait_b, tait_c;

        /**
         * The reference pressure
         */
        double reference_pressure;

        /**
         * The reference temperature
         */
        double reference_temperature;

        /**
         * The reference density
         */
        double reference_rho;

        /**
         * The reference isothermal bulk modulus
         */
        double reference_isothermal_bulk_modulus;

        /**
         * The reference K'
         * (first pressure derivative of the isothermal bulk modulus)
         */
        double reference_Kprime;

        /**
         * The reference K''
         * (second pressure derivative of the isothermal bulk modulus)
         */
        double Kdprime_0;

        /**
         * The reference thermal expansivity
         */
        double reference_thermal_expansivity;

        /**
         * The Einstein temperature
         */
        double einstein_temperature;

        /**
         * Parsed function that specifies the specific heat at P0 when using the Function
         * method.
         */
        Functions::ParsedFunction<1> reference_heat_capacity_function;

        /**
         * The constant viscosity
         */
        double eta;

        /**
         * The constant thermal conductivity.
         */
        double k_value;

    };

  }
}

#endif
