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
         * Specify which material model outputs the heating model requires
         * for computing the heating terms.
         */
        MaterialModel::MaterialProperties::Property
        get_required_properties () const override;

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
         *
         * This variable is read from the parameter file through the 'Limit
         * stress contribution to shear heating' parameter.
         */
        bool limit_stress;
        /**
         * This variable is read from the parameter file through a parameter called 'Cohesion for maximum shear stress'.
         */
        double cohesion;
        /**
         * This variable is read from the parameter file through a parameter called 'Friction angle for maximum shear stress'.
         */
        double friction_angle;
        MaterialModel::Rheology::DruckerPrager<dim> drucker_prager_plasticity;
    };


    /**
     * Additional output fields for the shear heating computation
     * to be added to the MaterialModel::MaterialModelOutputs structure
     * and filled in the MaterialModel::evaluate() function.
     *
     * This structure allows to modify the shear heating term by
     * multiplying it with a factor computed by the material model.
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

    /**
     * Additional output fields for the shear heating computation
     * to be added to the MaterialModel::MaterialModelOutputs structure
     * and filled in the MaterialModel::evaluate() function.
     *
     * This structure allows to prescribe the full shear heating term from
     * the material model.
     */
    template <int dim>
    class PrescribedShearHeatingOutputs : public MaterialModel::NamedAdditionalMaterialOutputs<dim>
    {
      public:
        PrescribedShearHeatingOutputs(const unsigned int n_points);

        std::vector<double> get_nth_output(const unsigned int idx) const override;

        /**
         * The viscous dissipation rate contributing to shear heating.
         * If this object is created and filled by the material model
         * it will replace the default viscous dissipation rate
         * computed by the shear heating model.
         */
        std::vector<double> prescribed_shear_heating_rates;
    };
  }
}


#endif
