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

#ifndef _aspect_material_model_composition_reaction_h
#define _aspect_material_model_composition_reaction_h

#include <aspect/material_model/interface.h>
#include <aspect/material_model/equation_of_state/linearized_incompressible.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * A material model that consists of globally constant values for all
     * material parameters except that the density decays linearly with the
     * temperature.
     *
     * The model is considered incompressible, following the definition
     * described in Interface::is_compressible. This is essentially the
     * material model used in the step-32 tutorial program.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class CompositionReaction : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
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

        /**
         * Create additional outputs
         */
        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const override;


      private:
        double reference_T;
        double eta;
        double composition_viscosity_prefactor_1;
        double composition_viscosity_prefactor_2;
        double thermal_viscosity_exponent;

        /**
         * The thermal conductivity.
         */
        double k_value;

        EquationOfState::LinearizedIncompressible<dim> equation_of_state;

        /**
         * Above this depth the compositional fields react: The first field
         * gets converted to the second field.
         */
        double reaction_depth;
    };

  }
}

#endif
