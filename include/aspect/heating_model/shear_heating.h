/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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


#ifndef _aspect_heating_model_shear_heating_h
#define _aspect_heating_model_shear_heating_h

#include <aspect/heating_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/rheology/drucker_prager.h>

namespace aspect
{
  namespace HeatingModel
  {
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
        void
        evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                  const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                  HeatingModel::HeatingModelOutputs &heating_model_outputs) const override;

        /**
         * Allow the heating model to attach additional material model outputs.
         */
        void
        create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &material_model_outputs) const override;

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

      private:
        /**
         * Parameters used for limiting the stress being used in the shear
         * heating computation. To prevent shear heating from becoming
         * unrealistically high, a Drucker-Prager yield criterion as given
         * in the DruckerPrager rheology model can be used to define a
         * maximum stress computed from the given cohesion and friction
         * angle.
         */
        bool limit_stress;
        double cohesion;
        double friction_angle;
        MaterialModel::Rheology::DruckerPrager<dim> drucker_prager_plasticity;
    };


    /**
     * Additional output fields for the shear heating computation
     * to be added to the MaterialModel::MaterialModelOutputs structure
     * and filled in the MaterialModel::evaluate() function.
     */
    template <int dim>
    class ShearHeatingOutputs : public MaterialModel::NamedAdditionalMaterialOutputs<dim>
    {
      public:
        ShearHeatingOutputs(const unsigned int n_points);

        std::vector<double> get_nth_output(const unsigned int idx) const override;

        /**
         * The fraction of the deformation work that is released as shear heating
         * rather than being converted into other forms of energy (such as, for
         * example, surface energy of grains). If it is set to 1, all deformation
         * work will go into shear heating.
         */
        std::vector<double> shear_heating_work_fractions;
    };
  }
}


#endif
