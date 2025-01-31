/*
  Copyright (C) 2020 - 2024 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_rheology_visco_plastic_h
#define _aspect_material_model_rheology_visco_plastic_h

#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <aspect/material_model/utilities.h>
#include <aspect/material_model/rheology/strain_dependent.h>
#include <aspect/material_model/rheology/friction_models.h>
#include <aspect/material_model/rheology/diffusion_creep.h>
#include <aspect/material_model/rheology/dislocation_creep.h>
#include <aspect/material_model/rheology/frank_kamenetskii.h>
#include <aspect/material_model/rheology/peierls_creep.h>
#include <aspect/material_model/rheology/constant_viscosity_prefactors.h>
#include <aspect/material_model/rheology/compositional_viscosity_prefactors.h>
#include <aspect/material_model/rheology/drucker_prager.h>
#include <aspect/material_model/rheology/elasticity.h>
#include <aspect/simulator_access.h>

#include<deal.II/fe/component_mask.h>

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * Additional output fields for the plastic parameters weakened (or hardened)
     * by strain to be added to the MaterialModel::MaterialModelOutputs structure
     * and filled in the MaterialModel::Interface::evaluate() function.
     */
    template <int dim>
    class PlasticAdditionalOutputs : public NamedAdditionalMaterialOutputs<dim>
    {
      public:
        PlasticAdditionalOutputs(const unsigned int n_points);

        std::vector<double> get_nth_output(const unsigned int idx) const override;

        /**
         * Cohesions at the evaluation points passed to
         * the instance of MaterialModel::Interface::evaluate() that fills
         * the current object.
         */
        std::vector<double> cohesions;

        /**
         * Internal angles of friction at the evaluation points passed to
         * the instance of MaterialModel::Interface::evaluate() that fills
         * the current object.
         */
        std::vector<double> friction_angles;

        /**
         * The plastic yield stress.
         */
        std::vector<double> yield_stresses;

        /**
         * The area where the viscous stress exceeds the plastic yield stress,
         * and viscosity is rescaled back to the yield envelope.
         */
        std::vector<double> yielding;

    };

    /**
     * A data structure with the output of calculate_isostrain_viscosities.
     */
    struct IsostrainViscosities
    {
      /**
       * The composition viscosity.
       */
      std::vector<double> composition_viscosities;

      /**
       * The composition yielding.
       */
      std::vector<bool> composition_yielding;

      /**
       * The current friction angle.
       */
      std::vector<double> current_friction_angles;

      /**
       * The current cohesion.
       */
      std::vector<double> current_cohesions;
    };

    namespace Rheology
    {

      template <int dim>
      class ViscoPlastic : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          ViscoPlastic();

          /**
           * This function calculates viscosities assuming that all the compositional fields
           * experience the same strain rate (isostrain).
           * If @p n_phase_transitions_per_composition points to a vector of
           * unsigned integers this is considered the number of phase transitions
           * for each compositional field and viscosity will be first computed on
           * each phase and then averaged for each compositional field.
           */
          IsostrainViscosities
          calculate_isostrain_viscosities ( const MaterialModel::MaterialModelInputs<dim> &in,
                                            const unsigned int i,
                                            const std::vector<double> &volume_fractions,
                                            const std::vector<double> &phase_function_values = std::vector<double>(),
                                            const std::vector<unsigned int> &n_phase_transitions_per_composition =
                                              std::vector<unsigned int>()) const;

          /**
           * A function that fills the viscosity derivatives in the
           * MaterialModelOutputs object that is handed over, if they exist.
           * Does nothing otherwise.
           * If @p n_phase_transitions_per_composition points to a vector of
           * unsigned integers this is considered the number of phase transitions
           * for each compositional field and viscosity will be first computed on
           * each phase and then averaged for each compositional field.
           */
          void compute_viscosity_derivatives(const unsigned int point_index,
                                             const std::vector<double> &volume_fractions,
                                             const std::vector<double> &composition_viscosities,
                                             const MaterialModel::MaterialModelInputs<dim> &in,
                                             MaterialModel::MaterialModelOutputs<dim> &out,
                                             const std::vector<double> &phase_function_values = std::vector<double>(),
                                             const std::vector<unsigned int> &n_phase_transitions_per_composition =
                                               std::vector<unsigned int>()) const;

          /**
           * A function that returns a ComponentMask that represents all compositional
           * fields that should be considered 'volumetric', that is representing a
           * physical proportion of the material, e.g. volume fraction of peridotite
           * (as opposed to non-volumetric quantities like the amount of finite-strain).
           */
          ComponentMask get_volumetric_composition_mask() const;

          /**
           * Declare the parameters this function takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           * If @p expected_n_phases_per_composition points to a vector of
           * unsigned integers this is considered the number of phases
           * for each compositional field and will be checked against the parsed
           * parameters.
           */
          void
          parse_parameters (ParameterHandler &prm,
                            const std::unique_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition = nullptr);

          /**
           * Create the additional material model outputs object that contains the
           * plastic outputs.
           */
          void
          create_plastic_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const;

          /**
           * A function that fills the plastic additional output in the
           * MaterialModelOutputs object that is handed over, if it exists.
           * Does nothing otherwise.
           */
          void fill_plastic_outputs(const unsigned int point_index,
                                    const std::vector<double> &volume_fractions,
                                    const bool plastic_yielding,
                                    const MaterialModel::MaterialModelInputs<dim> &in,
                                    MaterialModel::MaterialModelOutputs<dim> &out,
                                    const IsostrainViscosities &isostrain_viscosities) const;

          /**
           * Minimum strain rate used to stabilize the strain rate dependent rheology.
           */
          double min_strain_rate;

          /**
           * Enumeration for selecting which viscosity averaging scheme to use.
           */
          MaterialUtilities::CompositionalAveragingOperation viscosity_averaging;

          /**
           * Object for computing the strain dependence of the rheology model.
           */
          Rheology::StrainDependent<dim> strain_rheology;

          /**
           * Object for computing the friction dependence of the rheology model.
           */
          Rheology::FrictionModels<dim> friction_models;

          /**
           * Object for computing viscoelastic viscosities and stresses.
           */
          Rheology::Elasticity<dim> elastic_rheology;


        private:

          /**
           * Reference strain rate for the first non-linear iteration
           * in the first time step.
           */
          double ref_strain_rate;

          /**
           * Minimum and maximum viscosities used to improve the
           * stability of the rheology model.
           * These parameters contain one value per composition and phase (potentially the same value).
           */
          std::vector<double> minimum_viscosity;
          std::vector<double> maximum_viscosity;

          /**
           * Enumeration for selecting which type of viscous flow law to use.
           * Select between diffusion, dislocation, frank_kamenetskii or composite.
           */
          enum ViscosityScheme
          {
            diffusion,
            dislocation,
            frank_kamenetskii,
            composite
          } viscous_flow_law;

          /**
           * Enumeration for selecting which type of yield mechanism to use.
           * Select between Drucker Prager and stress limiter.
           */
          enum YieldScheme
          {
            stress_limiter,
            drucker_prager
          } yield_mechanism;

          /**
           * Whether to allow negative pressures to be used in the computation
           * of plastic yield stresses and viscosities. If false, the minimum
           * pressure in the plasticity formulation will be set to zero.
           */
          bool allow_negative_pressures_in_plasticity;

          /**
           * Whether to use the adiabatic pressure instead of the full pressure (default)
           * when calculating creep (diffusion, dislocation, and peierls) viscosity.
           * This may be helpful in models where the full pressure has an unusually
           * large negative value arising from large negative dynamic pressure,
           * resulting in solver convergence issue and in some cases a viscosity
           * of zero.
           */
          bool use_adiabatic_pressure_in_creep;

          /**
           * Whether to use the adiabatic pressure instead of the full pressure
           * when calculating the plastic yield stress.
           * This may be helpful in models where the full pressure has
           * large variations resulting in solver convergence issues.
           * Be aware that this setting will change the plastic shear band angle.
           */
          bool use_adiabatic_pressure_in_plasticity;

          /**
           * List of exponents controlling the behavior of the stress limiter
           * yielding mechanism.
           */
          std::vector<double> exponents_stress_limiter;

          /**
           * Temperature gradient added to temperature used in the flow law.
           */
          double adiabatic_temperature_gradient_for_viscosity;

          /**
           * Objects for computing viscous creep viscosities.
           */
          Rheology::DiffusionCreep<dim> diffusion_creep;
          Rheology::DislocationCreep<dim> dislocation_creep;
          std::unique_ptr<Rheology::FrankKamenetskii<dim>> frank_kamenetskii_rheology;

          /**
           * Whether to include Peierls creep in the constitutive formulation.
           */
          bool use_peierls_creep;

          /**
           * Object for computing Peierls creep viscosities.
           */
          std::unique_ptr<Rheology::PeierlsCreep<dim>> peierls_creep;

          /**
           * Object for computing the viscosity multiplied by a constant prefactor.
           * This multiplication step is done just prior to calculating the effective
           * viscoelastic viscosity or plastic viscosity.
           */
          Rheology::ConstantViscosityPrefactors<dim> constant_viscosity_prefactors;

          /**
           * Object for computing the viscosity multiplied by a given prefactor term.
           */
          Rheology::CompositionalViscosityPrefactors<dim> compositional_viscosity_prefactors;

          /*
           * Object for computing plastic stresses, viscosities, and additional outputs
           */
          Rheology::DruckerPrager<dim> drucker_prager_plasticity;

          /*
           * Input parameters for the drucker prager plasticity.
           */
          Rheology::DruckerPragerParameters drucker_prager_parameters;

      };
    }
  }
}
#endif
