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

#include <aspect/material_model/rheology/visco_plastic.h>
#include <aspect/material_model/utilities.h>
#include <aspect/utilities.h>
#include <aspect/newton.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace
    {
      std::vector<std::string> make_plastic_additional_outputs_names()
      {
        std::vector<std::string> names;
        names.emplace_back("current_cohesions");
        names.emplace_back("current_friction_angles");
        names.emplace_back("current_yield_stresses");
        names.emplace_back("plastic_yielding");
        return names;
      }
    }

    template <int dim>
    PlasticAdditionalOutputs<dim>::PlasticAdditionalOutputs(const unsigned int n_points)
      : NamedAdditionalMaterialOutputs<dim>(make_plastic_additional_outputs_names()),
        cohesions(n_points, numbers::signaling_nan<double>()),
        friction_angles(n_points, numbers::signaling_nan<double>()),
        yield_stresses(n_points, numbers::signaling_nan<double>()),
        yielding(n_points, numbers::signaling_nan<double>())
    {}



    template <int dim>
    std::vector<double>
    PlasticAdditionalOutputs<dim>::get_nth_output(const unsigned int idx) const
    {
      AssertIndexRange (idx, 4);
      switch (idx)
        {
          case 0:
            return cohesions;

          case 1:
            return friction_angles;

          case 2:
            return yield_stresses;

          case 3:
            return yielding;

          default:
            AssertThrow(false, ExcInternalError());
        }
      // We will never get here, so just return something
      return cohesions;
    }



    namespace Rheology
    {

      template <int dim>
      ViscoPlastic<dim>::ViscoPlastic ()
        = default;



      template <int dim>
      IsostrainViscosities
      ViscoPlastic<dim>::
      calculate_isostrain_viscosities (const MaterialModel::MaterialModelInputs<dim> &in,
                                       const unsigned int i,
                                       const std::vector<double> &volume_fractions,
                                       const std::vector<double> &phase_function_values,
                                       const std::vector<unsigned int> &n_phase_transitions_per_composition) const
      {
        IsostrainViscosities output_parameters;

        // Initialize or fill variables used to calculate viscosities
        output_parameters.composition_yielding.resize(volume_fractions.size(), false);
        output_parameters.composition_viscosities.resize(volume_fractions.size(), numbers::signaling_nan<double>());
        output_parameters.drucker_prager_parameters.resize(volume_fractions.size());
        output_parameters.dilation_lhs_terms.resize(volume_fractions.size(), numbers::signaling_nan<double>());
        output_parameters.dilation_rhs_terms.resize(volume_fractions.size(), numbers::signaling_nan<double>());

        // Assemble current and old stress tensor if elastic behavior is enabled
        SymmetricTensor<2, dim> stress_0_advected = numbers::signaling_nan<SymmetricTensor<2, dim>>();
        SymmetricTensor<2, dim> stress_old = numbers::signaling_nan<SymmetricTensor<2, dim>>();
        double elastic_shear_modulus = numbers::signaling_nan<double>();
        if (this->get_parameters().enable_elasticity)
          {
            // The first set of stresses in in.composition holds $\tau^{0adv}$
            // after the first advection nonlinear iteration,
            // which is when the viscosity matters for the Stokes system.
            const std::vector<unsigned int> &stress_field_indices = this->introspection().get_indices_for_fields_of_type(CompositionalFieldDescription::stress);
            stress_0_advected = (Utilities::Tensors::to_symmetric_tensor<dim>(&in.composition[i][stress_field_indices[0]],
                                                                              &in.composition[i][stress_field_indices[0]]+SymmetricTensor<2, dim>::n_independent_components));

            // The old stresses are only changed in the operator splitting step and have been advected into
            // the current timestep. They are stored in the second set of n_independent_components.
            stress_old = (Utilities::Tensors::to_symmetric_tensor<dim>(&in.composition[i][stress_field_indices[SymmetricTensor<2, dim>::n_independent_components]],
                                                                       &in.composition[i][stress_field_indices[SymmetricTensor<2, dim>::n_independent_components]]+SymmetricTensor<2, dim>::n_independent_components));


            // Average the compositional contributions to elastic_shear_moduli here and use
            // a volume-averaged shear modulus in the loop over the compositions below.
            // Otherwise it is implied that each material is acting independently
            // (different rotations, different stress changes), but this is inconsistent with storing only one stress tensor.
            // The averaging used is the same as for the viscosity.
            const std::vector<double> &elastic_shear_moduli = elastic_rheology.get_elastic_shear_moduli();
            elastic_shear_modulus = MaterialUtilities::average_value(volume_fractions, elastic_shear_moduli, viscosity_averaging);
          }

        // Use a specified "reference" strain rate if the strain rate is not yet available
        // during the very first nonlinear iteration or before. This is to avoid division by zero and
        // to have a better estimate of the resulting viscosity during this time.
        // During later iterations and timesteps the strain rate is capped by a minimum value min_strain_rate.
        const bool use_reference_strainrate = this->simulator_is_past_initialization() == false
                                              ||
                                              (this->get_timestep_number() == 0 &&
                                               this->get_nonlinear_iteration() == 0);

        double edot_ii;
        if (use_reference_strainrate)
          edot_ii = ref_strain_rate;
        else
          // Calculate the square root of the second moment invariant for the deviatoric strain rate tensor.
          edot_ii = std::max(std::sqrt(std::max(-second_invariant(deviator(in.strain_rate[i])), 0.)),
                             min_strain_rate);

        // Calculate viscosities for each of the individual compositional phases
        for (unsigned int j=0; j < volume_fractions.size(); ++j)
          {
            // Step 1: viscous behavior
            double non_yielding_viscosity = numbers::signaling_nan<double>();

            // Choice of activation volume depends on whether there is an adiabatic temperature
            // gradient used when calculating the viscosity. This allows the same activation volume
            // to be used in incompressible and compressible models.
            const double temperature_for_viscosity = (this->simulator_is_past_initialization())
                                                     ?
                                                     in.temperature[i] + adiabatic_temperature_gradient_for_viscosity*in.pressure[i]
                                                     :
                                                     this->get_adiabatic_conditions().temperature(in.position[i]);

            AssertThrow(temperature_for_viscosity != 0, ExcMessage(
                          "The temperature used in the calculation of the visco-plastic rheology is zero. "
                          "This is not allowed, because this value is used to divide through. It is probably "
                          "being caused by the temperature being zero somewhere in the model. The relevant "
                          "values for debugging are: temperature (" + Utilities::to_string(in.temperature[i]) +
                          "), adiabatic_temperature_gradient_for_viscosity ("
                          + Utilities::to_string(adiabatic_temperature_gradient_for_viscosity) + ") and pressure ("
                          + Utilities::to_string(in.pressure[i]) + ")."));
            {

              // Step 1a: compute viscosity from diffusion creep law, at least if it is going to be used

              // Determine whether to use the adiabatic pressure instead of the full pressure (default)
              // when calculating creep viscosity.
              double pressure_for_creep = in.pressure[i];

              if (use_adiabatic_pressure_in_creep)
                pressure_for_creep = this->get_adiabatic_conditions().pressure(in.position[i]);

              const double viscosity_diffusion
                = (viscous_flow_law != dislocation
                   ?
                   diffusion_creep.compute_viscosity(pressure_for_creep, temperature_for_viscosity, j,
                                                     phase_function_values,
                                                     n_phase_transitions_per_composition)
                   :
                   numbers::signaling_nan<double>());

              // Step 1b: compute viscosity from dislocation creep law
              const double viscosity_dislocation
                = (viscous_flow_law != diffusion
                   ?
                   dislocation_creep.compute_viscosity(edot_ii, pressure_for_creep, temperature_for_viscosity, j,
                                                       phase_function_values,
                                                       n_phase_transitions_per_composition)
                   :
                   numbers::signaling_nan<double>());

              // Step 1c: select which form of viscosity to use (diffusion, dislocation, their minimum or composite, or fk), and apply
              // pre-exponential weakening, if required.
              switch (viscous_flow_law)
                {
                  case diffusion:
                  {
                    non_yielding_viscosity = compositional_viscosity_prefactors.compute_viscosity(in, viscosity_diffusion, j, i, \
                                                                                                  CompositionalViscosityPrefactors<dim>::ModifiedFlowLaws::diffusion);
                    break;
                  }
                  case dislocation:
                  {
                    non_yielding_viscosity = compositional_viscosity_prefactors.compute_viscosity(in, viscosity_dislocation, j, i, \
                                                                                                  CompositionalViscosityPrefactors<dim>::ModifiedFlowLaws::dislocation);
                    break;
                  }
                  case frank_kamenetskii:
                  {
                    non_yielding_viscosity = frank_kamenetskii_rheology->compute_viscosity(in.temperature[i], j,
                                                                                           pressure_for_creep,
                                                                                           this->get_adiabatic_conditions().density(this->get_geometry_model().representative_point(0)),
                                                                                           this->get_gravity_model().gravity_vector(in.position[0]).norm());
                    break;
                  }
                  case composite:
                  {
                    const double scaled_viscosity_diffusion = compositional_viscosity_prefactors.compute_viscosity(in, viscosity_diffusion, j, i, \
                                                              CompositionalViscosityPrefactors<dim>::ModifiedFlowLaws::diffusion);
                    const double scaled_viscosity_dislocation = compositional_viscosity_prefactors.compute_viscosity(in, viscosity_dislocation, j, i, \
                                                                CompositionalViscosityPrefactors<dim>::ModifiedFlowLaws::dislocation);
                    non_yielding_viscosity = (scaled_viscosity_diffusion * scaled_viscosity_dislocation)/
                                             (scaled_viscosity_diffusion + scaled_viscosity_dislocation);
                    break;
                  }
                  case minimum_diffusion_dislocation:
                  {
                    const double scaled_viscosity_diffusion = compositional_viscosity_prefactors.compute_viscosity(in, viscosity_diffusion, j, i, \
                                                              CompositionalViscosityPrefactors<dim>::ModifiedFlowLaws::diffusion);
                    const double scaled_viscosity_dislocation = compositional_viscosity_prefactors.compute_viscosity(in, viscosity_dislocation, j, i, \
                                                                CompositionalViscosityPrefactors<dim>::ModifiedFlowLaws::dislocation);

                    non_yielding_viscosity = std::min(scaled_viscosity_diffusion, scaled_viscosity_dislocation);
                    break;
                  }
                  default:
                  {
                    AssertThrow(false, ExcNotImplemented());
                    break;
                  }
                }

              // Step 1d: compute the viscosity from the Peierls creep law and harmonically average with current viscosities
              if (use_peierls_creep)
                {
                  const double viscosity_peierls = peierls_creep->compute_viscosity(edot_ii, pressure_for_creep, temperature_for_viscosity, j,
                                                                                    phase_function_values,
                                                                                    n_phase_transitions_per_composition);
                  non_yielding_viscosity = (non_yielding_viscosity * viscosity_peierls) / (non_yielding_viscosity + viscosity_peierls);
                }
            }


            // Step 1e: multiply the viscosity by a constant (default value is 1)
            non_yielding_viscosity = constant_viscosity_prefactors.compute_viscosity(non_yielding_viscosity, j);

            // Step 2: calculate strain weakening factors for the cohesion, friction, and pre-yield viscosity
            // If no strain weakening is applied, the factors are 1.
            std::array<double, 3> weakening_factors = strain_rheology.compute_strain_weakening_factors(in.composition[i], j);

            if (strain_rheology.use_temperature_activated_strain_softening)
              weakening_factors = strain_rheology.apply_temperature_dependence_to_strain_weakening_factors(weakening_factors,temperature_for_viscosity,j);

            // Apply strain weakening to the viscous viscosity.
            non_yielding_viscosity *= weakening_factors[2];


            // Step 3: calculate the viscous stress magnitude
            // and strain rate. If requested compute visco-elastic contributions.
            double effective_edot_ii = edot_ii;

            if (this->get_parameters().enable_elasticity)
              {
                // Step 3a: calculate viscoelastic (effective) viscosity
                non_yielding_viscosity = elastic_rheology.calculate_viscoelastic_viscosity(non_yielding_viscosity, elastic_shear_modulus);

                if (use_reference_strainrate == false)
                  {
                    // Overwrite effective_edot_ii with a value that includes a term that accounts for
                    // elastic stress arising from a previous time step.
                    // If used, this variable is no longer the "true" strain rate, but is instead
                    // an effective value that enables the use of a standard isotropic viscosity
                    // formulation (i.e. where stress and strain are related by a scalar viscosity).
                    // The additional strain rate component is supported by a corresponding fictional body force.
                    Assert(std::isfinite(in.strain_rate[i].norm()),
                           ExcMessage("Invalid strain_rate in the MaterialModelInputs. This is likely because it was "
                                      "not filled by the caller."));

                    const SymmetricTensor<2, dim> effective_strain_rate = elastic_rheology.calculate_viscoelastic_strain_rate(in.strain_rate[i],
                                                                          stress_0_advected,
                                                                          stress_old,
                                                                          non_yielding_viscosity,
                                                                          elastic_shear_modulus);

                    effective_edot_ii = std::max(std::sqrt(std::max(-second_invariant(effective_strain_rate), 0.)),
                                                 min_strain_rate);
                  }
              }

            // Step 3b: calculate non yielding (viscous or viscous + elastic) stress magnitude
            const double non_yielding_stress = 2. * non_yielding_viscosity * effective_edot_ii;

            // Step 4a: calculate the strain-weakened friction and cohesion
            output_parameters.drucker_prager_parameters[j] = drucker_prager_plasticity.compute_drucker_prager_parameters(j,
                                                             phase_function_values,
                                                             n_phase_transitions_per_composition);
            output_parameters.drucker_prager_parameters[j].cohesion *= weakening_factors[0];
            output_parameters.drucker_prager_parameters[j].angle_internal_friction *= weakening_factors[1];

            // Step 4b: calculate the friction angle dependent on strain rate if specified
            // apply the strain rate dependence to the friction angle (including strain weakening if present)
            // Note: Maybe this should also be turned around to first apply strain rate dependence and then
            // the strain weakening to the dynamic friction angle. Didn't come up with a clear argument for
            // one order or the other.
            output_parameters.drucker_prager_parameters[j].angle_internal_friction = friction_models.compute_friction_angle(effective_edot_ii,
                                                                                     j,
                                                                                     output_parameters.drucker_prager_parameters[j].angle_internal_friction,
                                                                                     in.position[i]);

            // Step 5: plastic yielding

            // Determine if the pressure used in Drucker Prager plasticity will be capped at 0 (default).
            // This may be necessary in models without gravity and when the dynamic stresses are much higher
            // than the lithostatic pressure.

            double pressure_for_plasticity = in.pressure[i];

            if (use_adiabatic_pressure_in_plasticity)
              pressure_for_plasticity = this->get_adiabatic_conditions().pressure(in.position[i]);

            if (allow_negative_pressures_in_plasticity == false)
              pressure_for_plasticity = std::max(pressure_for_plasticity,0.0);

            // Step 5a: calculate the Drucker-Prager yield stress
            const double yield_stress = drucker_prager_plasticity.compute_yield_stress(pressure_for_plasticity,
                                                                                       output_parameters.drucker_prager_parameters[j]);

            // Step 5b: select if the yield viscosity is based on Drucker Prager or a stress limiter rheology
            double effective_viscosity = non_yielding_viscosity;
            switch (yield_mechanism)
              {
                case stress_limiter:
                {
                  //Step 5b-1: always rescale the viscosity back to the yield surface
                  const double viscosity_limiter = yield_stress / (2.0 * ref_strain_rate)
                                                   * std::pow((effective_edot_ii/ref_strain_rate),
                                                              1./exponents_stress_limiter[j] - 1.0);
                  effective_viscosity = 1. / ( 1./viscosity_limiter + 1./non_yielding_viscosity);
                  break;
                }
                case drucker_prager:
                {
                  // Step 5b-2: if the non-yielding stress is greater than the yield stress,
                  // rescale the viscosity back to yield surface
                  if (non_yielding_stress >= yield_stress)
                    {
                      // The following uses the effective_edot_ii
                      // (which has been modified for elastic effects, above),
                      // and calculates the effective viscosity over all active rheological elements
                      // assuming that the non-yielding viscosity is not strain rate dependent.
                      // TODO When dislocation creep or nonlinear plastic strain weakening
                      // is used, we need to iterate on the yield stress.
                      effective_viscosity = drucker_prager_plasticity.compute_viscosity(pressure_for_plasticity,
                                                                                        effective_edot_ii,
                                                                                        output_parameters.drucker_prager_parameters[j],
                                                                                        non_yielding_viscosity);
                      output_parameters.composition_yielding[j] = true;
                    }
                  break;
                }
                default:
                {
                  AssertThrow(false, ExcNotImplemented());
                  break;
                }
              }

            // Step 6: limit the viscosity with specified minimum and maximum bounds
            const double maximum_viscosity_for_composition = MaterialModel::MaterialUtilities::phase_average_value(
                                                               phase_function_values,
                                                               n_phase_transitions_per_composition,
                                                               maximum_viscosity,
                                                               j,
                                                               MaterialModel::MaterialUtilities::PhaseUtilities::logarithmic
                                                             );
            const double minimum_viscosity_for_composition = MaterialModel::MaterialUtilities::phase_average_value(
                                                               phase_function_values,
                                                               n_phase_transitions_per_composition,
                                                               minimum_viscosity,
                                                               j,
                                                               MaterialModel::MaterialUtilities::PhaseUtilities::logarithmic
                                                             );
            output_parameters.composition_viscosities[j] = std::min(std::max(effective_viscosity, minimum_viscosity_for_composition), maximum_viscosity_for_composition);

            // Compute the dilation terms if necessary.
            if (this->get_parameters().enable_prescribed_dilation == true)
              {
                if (output_parameters.composition_yielding[j] == true)
                  {
                    output_parameters.drucker_prager_parameters[j].angle_dilation *= weakening_factors[1];
                    const std::pair<double,double> dilation_terms = drucker_prager_plasticity.compute_dilation_terms_for_stokes_system (output_parameters.drucker_prager_parameters[j],
                                                                    non_yielding_viscosity,
                                                                    effective_edot_ii);
                    output_parameters.dilation_lhs_terms[j] = dilation_terms.first;
                    output_parameters.dilation_rhs_terms[j] = dilation_terms.second;
                  }
                else
                  {
                    output_parameters.dilation_lhs_terms[j] = 0;
                    output_parameters.dilation_rhs_terms[j] = 0;
                  }
              }
          }

        return output_parameters;
      }



      template <int dim>
      void
      ViscoPlastic<dim>::
      compute_viscosity_derivatives(const unsigned int i,
                                    const std::vector<double> &volume_fractions,
                                    const IsostrainViscosities &current_isostrain_values,
                                    const MaterialModel::MaterialModelInputs<dim> &in,
                                    MaterialModel::MaterialModelOutputs<dim> &out,
                                    const std::vector<double> &phase_function_values,
                                    const std::vector<unsigned int> &n_phase_transitions_per_composition) const
      {
        const std::shared_ptr<MaterialModel::MaterialModelDerivatives<dim>> derivatives
          = out.template get_additional_output_object<MaterialModel::MaterialModelDerivatives<dim>>();

        const bool enable_dilation = this->get_parameters().enable_prescribed_dilation;

        if (derivatives != nullptr)
          {
            // compute derivatives if necessary
            const unsigned int n_compositions = volume_fractions.size();

            std::vector<SymmetricTensor<2,dim>> composition_viscosity_derivatives_wrt_strain_rate(n_compositions);
            std::vector<double> composition_viscosity_derivatives_wrt_pressure(n_compositions);

            std::vector<SymmetricTensor<2,dim>> composition_dilation_derivatives_wrt_strain_rate(n_compositions);
            std::vector<double> composition_dilation_derivatives_wrt_pressure(n_compositions);

            const double finite_difference_accuracy = 1e-7;

            // A new material model inputs variable that uses the strain rate and pressure difference.
            MaterialModel::MaterialModelInputs<dim> in_derivatives = in;

            Assert(std::isfinite(in.strain_rate[i].norm()),
                   ExcMessage("Invalid strain_rate in the MaterialModelInputs. This is likely because it was "
                              "not filled by the caller."));

            const SymmetricTensor<2,dim> deviatoric_strain_rate = deviator(in.strain_rate[i]);

            // For each independent component, compute the derivative.
            for (unsigned int component = 0; component < SymmetricTensor<2,dim>::n_independent_components; ++component)
              {
                const TableIndices<2> strain_rate_indices = SymmetricTensor<2,dim>::unrolled_to_component_indices (component);

                const double strain_rate_difference = std::max(std::fabs(deviatoric_strain_rate[strain_rate_indices]), min_strain_rate)
                                                      * finite_difference_accuracy;
                const SymmetricTensor<2,dim> forward_strain_rate = deviatoric_strain_rate
                                                                   + strain_rate_difference * Utilities::nth_basis_for_symmetric_tensors<dim>(component);

                in_derivatives.strain_rate[i] = forward_strain_rate;

                const IsostrainViscosities forward_isostrain_values =
                  calculate_isostrain_viscosities(in_derivatives, i, volume_fractions,
                                                  phase_function_values, n_phase_transitions_per_composition);

                // For each composition of the independent component, compute the derivative.
                for (unsigned int composition_index = 0; composition_index < n_compositions; ++composition_index)
                  {
                    composition_viscosity_derivatives_wrt_strain_rate[composition_index][strain_rate_indices]
                      = ( forward_isostrain_values.composition_viscosities[composition_index] -
                          current_isostrain_values.composition_viscosities[composition_index]
                        ) / strain_rate_difference;

                    if (enable_dilation)
                      {
                        composition_dilation_derivatives_wrt_strain_rate[composition_index][strain_rate_indices]
                          = (
                              ( forward_isostrain_values.dilation_rhs_terms[composition_index] -
                                current_isostrain_values.dilation_rhs_terms[composition_index]
                              )
                              -
                              ( forward_isostrain_values.dilation_lhs_terms[composition_index] -
                                current_isostrain_values.dilation_lhs_terms[composition_index]
                              ) * in.pressure[i]
                            )
                            / strain_rate_difference;
                      }
                  }
              }

            // Now compute the derivative of the viscosity to the pressure
            const double pressure_difference = std::fabs(in.pressure[i]) * finite_difference_accuracy;
            const double forward_pressure = in.pressure[i] + pressure_difference;

            in_derivatives.pressure[i] = forward_pressure;

            // Modify the in_derivatives object again to take the original strain rate.
            in_derivatives.strain_rate[i] = in.strain_rate[i];

            const IsostrainViscosities forward_isostrain_values =
              calculate_isostrain_viscosities(in_derivatives, i, volume_fractions,
                                              phase_function_values, n_phase_transitions_per_composition);

            for (unsigned int composition_index = 0; composition_index < n_compositions; ++composition_index)
              {
                composition_viscosity_derivatives_wrt_pressure[composition_index]
                  = pressure_difference > 0.0 ?
                    ( forward_isostrain_values.composition_viscosities[composition_index] -
                      current_isostrain_values.composition_viscosities[composition_index]
                    ) / pressure_difference :
                    0.0;

                if (enable_dilation)
                  {
                    composition_dilation_derivatives_wrt_pressure[composition_index]
                      = pressure_difference > 0.0
                        ?
                        (
                          ( forward_isostrain_values.dilation_rhs_terms[composition_index] -
                            current_isostrain_values.dilation_rhs_terms[composition_index]
                          )
                          -
                          ( forward_isostrain_values.dilation_lhs_terms[composition_index] -
                            current_isostrain_values.dilation_lhs_terms[composition_index]
                          ) * in.pressure[i]
                        )
                        / pressure_difference
                        :
                        0.0;
                  }
              }

            double viscosity_averaging_p = 0; // Geometric
            if (viscosity_averaging == MaterialUtilities::harmonic)
              viscosity_averaging_p = -1;
            if (viscosity_averaging == MaterialUtilities::arithmetic)
              viscosity_averaging_p = 1;
            if (viscosity_averaging == MaterialUtilities::maximum_composition)
              viscosity_averaging_p = 1000;

            derivatives->viscosity_derivative_wrt_strain_rate[i] =
              Utilities::derivative_of_weighted_p_norm_average(out.viscosities[i],
                                                               volume_fractions,
                                                               current_isostrain_values.composition_viscosities,
                                                               composition_viscosity_derivatives_wrt_strain_rate,
                                                               viscosity_averaging_p);
            derivatives->viscosity_derivative_wrt_pressure[i] =
              Utilities::derivative_of_weighted_p_norm_average(out.viscosities[i],
                                                               volume_fractions,
                                                               current_isostrain_values.composition_viscosities,
                                                               composition_viscosity_derivatives_wrt_pressure,
                                                               viscosity_averaging_p);

            if (enable_dilation)
              {
                std::vector<double> dilation_values(n_compositions);
                for (unsigned int composition_index = 0; composition_index < n_compositions; ++composition_index)
                  dilation_values[composition_index] = current_isostrain_values.dilation_rhs_terms[composition_index] -
                                                       current_isostrain_values.dilation_lhs_terms[composition_index] * in.pressure[i];

                double dummy_averaged_parameter = numbers::signaling_nan<double>();

                // always use arithmetic average for plastic dilation terms
                derivatives->dilation_derivative_wrt_strain_rate[i] =
                  Utilities::derivative_of_weighted_p_norm_average(dummy_averaged_parameter,
                                                                   volume_fractions,
                                                                   dilation_values,
                                                                   composition_dilation_derivatives_wrt_strain_rate,
                                                                   /*arithmetic average =*/ 1);

                derivatives->dilation_derivative_wrt_pressure[i] =
                  Utilities::derivative_of_weighted_p_norm_average(dummy_averaged_parameter,
                                                                   volume_fractions,
                                                                   dilation_values,
                                                                   composition_dilation_derivatives_wrt_pressure,
                                                                   /*arithmetic average =*/ 1);
              }
          }
      }



      template <int dim>
      ComponentMask
      ViscoPlastic<dim>::
      get_volumetric_composition_mask() const
      {
        // Store which components to exclude during the volume fraction computation.
        ComponentMask composition_mask = strain_rheology.get_strain_composition_mask();

        if (this->get_parameters().enable_elasticity)
          {
            // Set a mask (false) for the 2*n_independent_components fields representing the viscoelastic
            // stress tensor components so that they are excluded from the volume fraction computation.
            const std::vector<unsigned int> &stress_field_indices = this->introspection().get_indices_for_fields_of_type(CompositionalFieldDescription::stress);
            for (unsigned int field_index: stress_field_indices)
              composition_mask.set(field_index, false);
          }

        return composition_mask;
      }



      template <int dim>
      void
      ViscoPlastic<dim>::declare_parameters (ParameterHandler &prm)
      {
        Rheology::StrainDependent<dim>::declare_parameters (prm);

        Rheology::FrictionModels<dim>::declare_parameters (prm);

        Rheology::Elasticity<dim>::declare_parameters (prm);

        // Reference and minimum/maximum values
        prm.declare_entry ("Minimum strain rate", "1.0e-20", Patterns::Double (0.),
                           "Stabilizes strain dependent viscosity. Units: \\si{\\per\\second}.");
        prm.declare_entry ("Reference strain rate","1.0e-15",Patterns::Double (0.),
                           "Reference strain rate for first time step. Units: \\si{\\per\\second}.");
        prm.declare_entry ("Minimum viscosity", "1e17", Patterns::Anything(),
                           "Lower cutoff for effective viscosity. Units: $\\text{Pa}\\text{s}$. "
                           "List with as many components as active "
                           "compositional fields (material data is assumed to "
                           "be in order with the ordering of the fields). ");
        prm.declare_entry ("Maximum viscosity", "1e28", Patterns::Anything(),
                           "Upper cutoff for effective viscosity. Units: $\\text{Pa}\\text{s}$. "
                           "List with as many components as active "
                           "compositional fields (material data is assumed to "
                           "be in order with the ordering of the fields). ");

        // Rheological parameters
        prm.declare_entry ("Viscosity averaging scheme", "harmonic",
                           Patterns::Selection("arithmetic|harmonic|geometric|maximum composition"),
                           "When more than one compositional field is present at a point "
                           "with different viscosities, we need to come up with an average "
                           "viscosity at that point.  Select a weighted harmonic, arithmetic, "
                           "geometric, or maximum composition.");
        prm.declare_entry ("Viscous flow law", "composite",
                           Patterns::Selection("diffusion|dislocation|frank kamenetskii|composite|minimum diffusion dislocation"),
                           "Select what type of viscosity law to use between the options diffusion, "
                           "dislocation, frank kamenetskii, composite and the minimum of diffusion and dislocation. "
                           "Soon there will be an option "
                           "to select a specific flow law for each assigned composition, "
                           "When the full strain rate is used to compute each viscosity "
                           "instead of the properly partitioned diffusion and dislocation creep "
                           "strain rates, it is recommended to take the minimum of the creep "
                           "viscosities instead of their composite.");
        prm.declare_entry ("Yield mechanism", "drucker",
                           Patterns::Selection("drucker|limiter"),
                           "Select what type of yield mechanism to use between Drucker Prager "
                           "and stress limiter options.");
        prm.declare_entry ("Allow negative pressures in plasticity", "false",
                           Patterns::Bool (),
                           "Whether to allow negative pressures to be used in the computation "
                           "of plastic yield stresses and viscosities. Setting this parameter "
                           "to true may be advantageous in models without gravity where the "
                           "dynamic stresses are much higher than the lithostatic pressure. "
                           "If false, the minimum pressure in the plasticity formulation will "
                           "be set to zero.");
        prm.declare_entry ("Use adiabatic pressure in creep viscosity", "false",
                           Patterns::Bool (),
                           "Whether to use the adiabatic pressure instead of the full "
                           "pressure (default) when calculating viscous creep. "
                           "This may be helpful in models where the "
                           "full pressure has an unusually large negative value arising from "
                           "large negative dynamic pressure, resulting in solver convergence "
                           "issue and in some cases a viscosity of zero.");
        prm.declare_entry ("Use adiabatic pressure in plasticity", "false",
                           Patterns::Bool (),
                           "Whether to use the adiabatic pressure instead of the full "
                           "pressure when calculating plastic yield stress. "
                           "This may be helpful in models where the "
                           "full pressure has unusually large variations, resulting "
                           "in solver convergence issues. Be aware that this setting "
                           "will change the plastic shear band angle.");

        // Diffusion creep parameters
        Rheology::DiffusionCreep<dim>::declare_parameters(prm);

        // Dislocation creep parameters
        Rheology::DislocationCreep<dim>::declare_parameters(prm);

        // Frank-Kamenetskii viscosity parameters
        Rheology::FrankKamenetskii<dim>::declare_parameters(prm);

        // Peierls creep parameters
        Rheology::PeierlsCreep<dim>::declare_parameters(prm);

        prm.declare_entry ("Include Peierls creep", "false",
                           Patterns::Bool (),
                           "Whether to include Peierls creep in the rheological formulation.");

        // Constant viscosity prefactor parameters
        Rheology::ConstantViscosityPrefactors<dim>::declare_parameters(prm);

        // Variable viscosity prefactor parameters
        Rheology::CompositionalViscosityPrefactors<dim>::declare_parameters(prm);

        // Drucker Prager plasticity parameters
        Rheology::DruckerPrager<dim>::declare_parameters(prm);

        // Stress limiter parameters
        prm.declare_entry ("Stress limiter exponents", "1.0",
                           Patterns::List(Patterns::Double (0.)),
                           "List of stress limiter exponents, $n_{\\text{lim}}$, "
                           "for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. Units: none.");

        // Temperature gradient in viscosity laws to include an adiabat (note units of K/Pa)
        prm.declare_entry ("Adiabat temperature gradient for viscosity", "0.0", Patterns::Double (0.),
                           "Add an adiabatic temperature gradient to the temperature used in the flow law "
                           "so that the activation volume is consistent with what one would use in a "
                           "earth-like (compressible) model. Default is set so this is off. "
                           "Note that this is a linear approximation of the real adiabatic gradient, which "
                           "is okay for the upper mantle, but is not really accurate for the lower mantle. "
                           "Using a pressure gradient of 32436 Pa/m, then a value of "
                           "0.3 K/km = 0.0003 K/m = 9.24e-09 K/Pa gives an earth-like adiabat."
                           "Units: \\si{\\kelvin\\per\\pascal}.");
      }



      template <int dim>
      void
      ViscoPlastic<dim>::parse_parameters (ParameterHandler &prm,
                                           const std::unique_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition)
      {
        strain_rheology.initialize_simulator (this->get_simulator());
        strain_rheology.parse_parameters(prm);

        friction_models.initialize_simulator (this->get_simulator());
        friction_models.parse_parameters(prm);

        if (this->get_parameters().enable_elasticity)
          {
            elastic_rheology.initialize_simulator (this->get_simulator());
            elastic_rheology.parse_parameters(prm);
          }

        // Reference and minimum/maximum values
        min_strain_rate = prm.get_double("Minimum strain rate");
        ref_strain_rate = prm.get_double("Reference strain rate");
        AssertThrow(ref_strain_rate >= min_strain_rate,
                    ExcMessage("The reference strain rate for the viscoplastic material model should be larger than the minimum strain rate."));

        // Retrieve the list of composition names
        std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();

        // Retrieve the list of names of fields that represent chemical compositions, and not, e.g.,
        // plastic strain
        std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();

        // Establish that a background field is required here
        compositional_field_names.insert(compositional_field_names.begin(), "background");
        chemical_field_names.insert(chemical_field_names.begin(), "background");

        Utilities::MapParsing::Options options(chemical_field_names, "Minimum viscosity");
        options.list_of_allowed_keys = compositional_field_names;
        options.allow_multiple_values_per_key = true;
        if (expected_n_phases_per_composition)
          {
            options.n_values_per_key = *expected_n_phases_per_composition;

            // check_values_per_key is required to be true to duplicate single values
            // if they are to be used for all phases associated with a given key.
            options.check_values_per_key = true;
          }

        minimum_viscosity = Utilities::MapParsing::parse_map_to_double_array (prm.get("Minimum viscosity"),
                                                                              options);

        options.property_name = "Maximum viscosity";
        maximum_viscosity = Utilities::MapParsing::parse_map_to_double_array (prm.get("Maximum viscosity"),
                                                                              options);

        Assert(maximum_viscosity.size() == minimum_viscosity.size(),
               ExcMessage("The input parameters 'Maximum viscosity' and 'Minimum viscosity' should have the same number of entries."));
        for (auto p1 = maximum_viscosity.begin(), p2 = minimum_viscosity.begin();
             p1 != maximum_viscosity.end(); ++p1, ++p2)
          AssertThrow(*p1 >= *p2, ExcMessage("Maximum viscosity should be larger or equal to the minimum viscosity."));

        viscosity_averaging = MaterialUtilities::parse_compositional_averaging_operation ("Viscosity averaging scheme",
                              prm);

        // Rheological parameters
        if (prm.get ("Viscous flow law") == "composite")
          viscous_flow_law = composite;
        else if (prm.get ("Viscous flow law") == "diffusion")
          viscous_flow_law = diffusion;
        else if (prm.get ("Viscous flow law") == "dislocation")
          viscous_flow_law = dislocation;
        else if (prm.get ("Viscous flow law") == "frank kamenetskii")
          viscous_flow_law = frank_kamenetskii;
        else if (prm.get ("Viscous flow law") == "minimum diffusion dislocation")
          viscous_flow_law = minimum_diffusion_dislocation;
        else
          AssertThrow(false, ExcMessage("Not a valid viscous flow law"));

        if (prm.get("Yield mechanism") == "drucker")
          yield_mechanism = drucker_prager;
        else if (prm.get ("Yield mechanism") == "limiter")
          yield_mechanism = stress_limiter;
        else
          AssertThrow(false, ExcMessage("Not a valid yield mechanism."));

        AssertThrow(this->get_parameters().enable_elasticity == false || yield_mechanism == drucker_prager,
                    ExcMessage("Elastic behavior is only tested with the "
                               "'drucker prager' plasticity option."));

        allow_negative_pressures_in_plasticity = prm.get_bool ("Allow negative pressures in plasticity");
        use_adiabatic_pressure_in_plasticity = prm.get_bool("Use adiabatic pressure in plasticity");
        use_adiabatic_pressure_in_creep = prm.get_bool("Use adiabatic pressure in creep viscosity");

        if (this->get_parameters().enable_prescribed_dilation)
          AssertThrow(allow_negative_pressures_in_plasticity == true &&
                      use_adiabatic_pressure_in_plasticity == false,
                      ExcMessage("Currently, plastic dilation requires that the plastic rheology "
                                 "depends on the true value of pressure."));

        // Diffusion creep parameters
        diffusion_creep.initialize_simulator (this->get_simulator());
        diffusion_creep.parse_parameters(prm, expected_n_phases_per_composition);

        // Dislocation creep parameters
        dislocation_creep.initialize_simulator (this->get_simulator());
        dislocation_creep.parse_parameters(prm, expected_n_phases_per_composition);

        // Frank Kamenetskii viscosity parameters
        if (viscous_flow_law == frank_kamenetskii)
          {
            frank_kamenetskii_rheology = std::make_unique<Rheology::FrankKamenetskii<dim>>();
            frank_kamenetskii_rheology->initialize_simulator (this->get_simulator());
            frank_kamenetskii_rheology->parse_parameters(prm);
          }

        // Peierls creep parameters
        use_peierls_creep = prm.get_bool ("Include Peierls creep");
        if (use_peierls_creep)
          {
            peierls_creep = std::make_unique<Rheology::PeierlsCreep<dim>>();
            peierls_creep->initialize_simulator (this->get_simulator());
            peierls_creep->parse_parameters(prm, expected_n_phases_per_composition);
          }

        // Constant viscosity prefactor parameters
        constant_viscosity_prefactors.initialize_simulator (this->get_simulator());
        constant_viscosity_prefactors.parse_parameters(prm);

        compositional_viscosity_prefactors.initialize_simulator (this->get_simulator());
        compositional_viscosity_prefactors.parse_parameters(prm);

        // Plasticity parameters
        drucker_prager_plasticity.initialize_simulator (this->get_simulator());
        drucker_prager_plasticity.parse_parameters(prm, expected_n_phases_per_composition);

        // Stress limiter parameter (does not allow for phases per composition)
        options.property_name = "Stress limiter exponents";
        options.allow_multiple_values_per_key = false;
        options.check_values_per_key = false;

        exponents_stress_limiter = Utilities::MapParsing::parse_map_to_double_array (prm.get("Stress limiter exponents"),
                                                                                     options);

        // Include an adiabat temperature gradient in flow laws
        adiabatic_temperature_gradient_for_viscosity = prm.get_double("Adiabat temperature gradient for viscosity");
        if (this->get_heating_model_manager().adiabatic_heating_enabled())
          AssertThrow (adiabatic_temperature_gradient_for_viscosity == 0.0,
                       ExcMessage("If adiabatic heating is enabled you should not add another adiabatic gradient"
                                  "to the temperature for computing the viscosity, because the ambient"
                                  "temperature profile already includes the adiabatic gradient."));

      }



      template <int dim>
      void
      ViscoPlastic<dim>::create_plastic_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        if (out.template has_additional_output_object<PlasticAdditionalOutputs<dim>>() == false)
          {
            const unsigned int n_points = out.n_evaluation_points();
            out.additional_outputs.push_back(
              std::make_unique<PlasticAdditionalOutputs<dim>> (n_points));
          }
      }

      template <int dim>
      void
      ViscoPlastic<dim>::
      fill_plastic_outputs(const unsigned int i,
                           const std::vector<double> &volume_fractions,
                           const bool plastic_yielding,
                           const MaterialModel::MaterialModelInputs<dim> &in,
                           MaterialModel::MaterialModelOutputs<dim> &out,
                           const IsostrainViscosities &isostrain_viscosities) const
      {
        const std::shared_ptr<PlasticAdditionalOutputs<dim>> plastic_out
          = out.template get_additional_output_object<PlasticAdditionalOutputs<dim>>();

        if (plastic_out != nullptr)
          {
            AssertThrow(!std::isnan(out.viscosities[0]),
                        ExcMessage("The PlasticAdditionalOutputs cannot be filled when the viscosity has not been computed."));

            plastic_out->cohesions[i] = 0;
            plastic_out->friction_angles[i] = 0;
            plastic_out->yield_stresses[i] = 0;
            plastic_out->yielding[i] = plastic_yielding ? 1 : 0;

            double pressure_for_plasticity = in.pressure[i];

            if (use_adiabatic_pressure_in_plasticity)
              pressure_for_plasticity = this->get_adiabatic_conditions().pressure(in.position[i]);

            if (allow_negative_pressures_in_plasticity == false)
              pressure_for_plasticity = std::max(pressure_for_plasticity, 0.0);

            // average over the volume volume fractions
            for (unsigned int j = 0; j < volume_fractions.size(); ++j)
              {
                const Rheology::DruckerPragerParameters &drucker_prager_parameters = isostrain_viscosities.drucker_prager_parameters[j];

                plastic_out->cohesions[i] += volume_fractions[j] * drucker_prager_parameters.cohesion;
                // Also convert radians to degrees
                plastic_out->friction_angles[i] += constants::radians_to_degree * volume_fractions[j] * drucker_prager_parameters.angle_internal_friction;

                plastic_out->yield_stresses[i] += volume_fractions[j] * drucker_prager_plasticity.compute_yield_stress(pressure_for_plasticity,
                                                  drucker_prager_parameters);
              }
          }
      }

    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  template class PlasticAdditionalOutputs<dim>; \
  \
  namespace Rheology \
  { \
    template class ViscoPlastic<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
