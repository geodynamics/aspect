/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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

#include <aspect/material_model/visco_plastic.h>
#include <aspect/utilities.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/signaling_nan.h>
#include <aspect/newton.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>

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
        names.emplace_back("plastic_yielding");
        return names;
      }
    }

    template <int dim>
    bool
    ViscoPlastic<dim>::
    is_yielding (const double &pressure,
                 const double &temperature,
                 const std::vector<double> &composition,
                 const SymmetricTensor<2,dim> &strain_rate) const
    {
      /* The following returns whether or not the material is plastically yielding
       * as documented in evaluate.
       */
      bool plastic_yielding = false;

      MaterialModel::MaterialModelInputs <dim> in (in.n_evaluation_points(), this->n_compositional_fields());
      unsigned int i = 0;

      in.pressure[i] = pressure;
      in.temperature[i] = temperature;
      in.composition[i] = composition;
      in.strain_rate[i] = strain_rate;

      const std::vector<double> volume_fractions = MaterialUtilities::compute_volume_fractions(composition, get_volumetric_composition_mask());

      const std::pair<std::vector<double>, std::vector<bool> > calculate_viscosities =
        calculate_isostrain_viscosities(in, i, volume_fractions, viscous_flow_law, yield_mechanism);

      std::vector<double>::const_iterator max_composition = std::max_element(volume_fractions.begin(),volume_fractions.end());
      plastic_yielding = calculate_viscosities.second[std::distance(volume_fractions.begin(),max_composition)];

      return plastic_yielding;
    }



    template <int dim>
    bool
    ViscoPlastic<dim>::
    is_yielding(const MaterialModelInputs<dim> &in) const
    {
      Assert(in.n_evaluation_points() == 1, ExcInternalError());

      const std::vector<double> volume_fractions = MaterialUtilities::compute_volume_fractions(in.composition[0], get_volumetric_composition_mask());

      /* The following handles phases in a similar way as in the 'evaluate' function.
      * Results then enter the calculation of plastic yielding.
      */
      std::vector<double> phase_function_values(phase_function.n_phase_transitions(), 0.0);

      if (phase_function.n_phase_transitions() > 0)
        {
          const double gravity_norm = this->get_gravity_model().gravity_vector(in.position[0]).norm();

          double reference_density;
          if (this->get_adiabatic_conditions().is_initialized())
            {
              reference_density = this->get_adiabatic_conditions().density(in.position[0]);
            }
          else
            {
              EquationOfStateOutputs<dim> eos_outputs_all_phases (this->n_compositional_fields()+1+phase_function.n_phase_transitions());
              equation_of_state.evaluate(in, 0, eos_outputs_all_phases);
              reference_density = eos_outputs_all_phases.densities[0];
            }

          MaterialUtilities::PhaseFunctionInputs<dim> phase_inputs(in.temperature[0],
                                                                   in.pressure[0],
                                                                   this->get_geometry_model().depth(in.position[0]),
                                                                   gravity_norm*reference_density,
                                                                   numbers::invalid_unsigned_int);


          for (unsigned int j=0; j < phase_function.n_phase_transitions(); j++)
            {
              phase_inputs.phase_index = j;
              phase_function_values[j] = phase_function.compute_value(phase_inputs);
            }
        }

      /* The following returns whether or not the material is plastically yielding
       * as documented in evaluate.
       */

      const std::pair<std::vector<double>, std::vector<bool>> calculate_viscosities =
                                                             calculate_isostrain_viscosities(in, 0, volume_fractions, viscous_flow_law,
                                                                 yield_mechanism, phase_function_values);

      std::vector<double>::const_iterator max_composition = std::max_element(volume_fractions.begin(), volume_fractions.end());
      const bool plastic_yielding = calculate_viscosities.second[std::distance(volume_fractions.begin(), max_composition)];

      return plastic_yielding;
    }



    template <int dim>
    PlasticAdditionalOutputs<dim>::PlasticAdditionalOutputs (const unsigned int n_points)
      :
      NamedAdditionalMaterialOutputs<dim>(make_plastic_additional_outputs_names()),
      cohesions(n_points, numbers::signaling_nan<double>()),
      friction_angles(n_points, numbers::signaling_nan<double>()),
      yielding(n_points, numbers::signaling_nan<double>())
    {}

    template <int dim>
    std::vector<double>
    PlasticAdditionalOutputs<dim>::get_nth_output(const unsigned int idx) const
    {
      AssertIndexRange (idx, 3);
      switch (idx)
        {
          case 0:
            return cohesions;

          case 1:
            return friction_angles;

          case 2:
            return yielding;

          default:
            AssertThrow(false, ExcInternalError());
        }
      // We will never get here, so just return something
      return cohesions;
    }



    template <int dim>
    std::pair<std::vector<double>, std::vector<bool> >
    ViscoPlastic<dim>::
    calculate_isostrain_viscosities (const MaterialModel::MaterialModelInputs<dim> &in,
                                     const unsigned int i,
                                     const std::vector<double> &volume_fractions,
                                     const ViscosityScheme &viscous_type,
                                     const YieldScheme &yield_type,
                                     const std::vector<double> &phase_function_values) const
    {
      // Initialize or fill variables used to calculate viscosities
      std::vector<bool> composition_yielding(volume_fractions.size(), false);
      std::vector<double> composition_viscosities(volume_fractions.size(), numbers::signaling_nan<double>());

      // The first time this function is called (first iteration of first time step)
      // a specified "reference" strain rate is used as the returned value would
      // otherwise be zero.
      const bool use_reference_strainrate = (this->get_timestep_number() == 0) &&
                                            (in.strain_rate[i].norm() <= std::numeric_limits<double>::min());

      double edot_ii;
      if (use_reference_strainrate)
        edot_ii = ref_strain_rate;
      else
        // Calculate the square root of the second moment invariant for the deviatoric strain rate tensor.
        edot_ii = std::max(std::sqrt(std::fabs(second_invariant(deviator(in.strain_rate[i])))),
                           min_strain_rate);

      // Calculate viscosities for each of the individual compositional phases
      for (unsigned int j=0; j < volume_fractions.size(); ++j)
        {
          // Step 1: viscous behavior

          // Choice of activation volume depends on whether there is an adiabatic temperature
          // gradient used when calculating the viscosity. This allows the same activation volume
          // to be used in incompressible and compressible models.
          const double temperature_for_viscosity = in.temperature[i] + adiabatic_temperature_gradient_for_viscosity*in.pressure[i];
          AssertThrow(temperature_for_viscosity != 0, ExcMessage(
                        "The temperature used in the calculation of the visco-plastic rheology is zero. "
                        "This is not allowed, because this value is used to divide through. It is probably "
                        "being caused by the temperature being zero somewhere in the model. The relevant "
                        "values for debugging are: temperature (" + Utilities::to_string(in.temperature[i]) +
                        "), adiabatic_temperature_gradient_for_viscosity ("
                        + Utilities::to_string(adiabatic_temperature_gradient_for_viscosity) + ") and pressure ("
                        + Utilities::to_string(in.pressure[i]) + ")."));

          // Step 1a: compute viscosity from diffusion creep law
          const double viscosity_diffusion = diffusion_creep.compute_viscosity(in.pressure[i], temperature_for_viscosity, j,
                                                                               phase_function_values,
                                                                               phase_function.n_phase_transitions_for_each_composition());

          // Step 1b: compute viscosity from dislocation creep law
          const double viscosity_dislocation = dislocation_creep.compute_viscosity(edot_ii, in.pressure[i], temperature_for_viscosity, j,
                                                                                   phase_function_values,
                                                                                   phase_function.n_phase_transitions_for_each_composition());

          // Step 1c: select what form of viscosity to use (diffusion, dislocation, fk, or composite)
          double viscosity_pre_yield = 0.0;
          switch (viscous_type)
            {
              case diffusion:
              {
                viscosity_pre_yield = viscosity_diffusion;
                break;
              }
              case dislocation:
              {
                viscosity_pre_yield = viscosity_dislocation;
                break;
              }
              case frank_kamenetskii:
              {
                viscosity_pre_yield = frank_kamenetskii_rheology->compute_viscosity(in.temperature[i], j);
                break;
              }
              case composite:
              {
                viscosity_pre_yield = (viscosity_diffusion * viscosity_dislocation)/
                                      (viscosity_diffusion + viscosity_dislocation);
                break;
              }
              default:
              {
                AssertThrow(false, ExcNotImplemented());
                break;
              }
            }

          // Step 1d: compute viscosity from Peierls creep law and harmonically average with current viscosities
          if (use_peierls_creep)
            {
              const double viscosity_peierls = peierls_creep->compute_viscosity(edot_ii, in.pressure[i], temperature_for_viscosity, j);
              viscosity_pre_yield = (viscosity_pre_yield * viscosity_peierls) / (viscosity_pre_yield + viscosity_peierls);
            }

          // Step 1e: multiply the viscosity by a constant (default value is 1)
          viscosity_pre_yield = constant_viscosity_prefactors.compute_viscosity(viscosity_pre_yield, j);

          // Step 2: calculate the viscous stress magnitude
          // and strain rate. If requested compute visco-elastic contributions.
          const double dte = elastic_rheology.elastic_timestep();
          const std::vector<double> &elastic_shear_moduli = elastic_rheology.get_elastic_shear_moduli();
          const double current_edot_ii = Utilities::compute_current_edot_ii(j,
                                                                            in.composition[i],
                                                                            ref_strain_rate,
                                                                            min_strain_rate,
                                                                            in.strain_rate[i],
                                                                            elastic_shear_moduli,
                                                                            use_elasticity,
                                                                            use_reference_strainrate,
                                                                            dte);


          double current_stress = numbers::signaling_nan<double>();
          if (use_elasticity == false)
            current_stress = 2. * viscosity_pre_yield * current_edot_ii;
          else
            {
              if (use_reference_strainrate == false)
                {
                  // Step 2a: calculate viscoelastic (effective) viscosity
                  viscosity_pre_yield = elastic_rheology.calculate_viscoelastic_viscosity(viscosity_pre_yield,
                                                                                          elastic_shear_moduli[j]);

                  // Step 2b: calculate current (viscous + elastic) stress magnitude
                  current_stress = viscosity_pre_yield * current_edot_ii;
                }
            }

          // Step 3: strain weakening

          // Step 3a: calculate strain weakening factors for the cohesion, friction, and pre-yield viscosity
          // If no brittle and/or viscous strain weakening is applied, the factors are 1.
          const std::array<double, 3> weakening_factors = strain_rheology.compute_strain_weakening_factors(j, in.composition[i]);

          // Step 3b: calculate weakened friction, cohesion, and pre-yield viscosity
          const double current_cohesion = drucker_prager_parameters.cohesions[j] * weakening_factors[0];
          const double current_friction = drucker_prager_parameters.angles_internal_friction[j] * weakening_factors[1];
          viscosity_pre_yield *= weakening_factors[2];

          // Step 4: plastic yielding

          // Step 4a: calculate Drucker-Prager yield stress
          const double yield_stress = drucker_prager_plasticity.compute_yield_stress(current_cohesion,
                                                                                     current_friction,
                                                                                     std::max(in.pressure[i], 0.0),
                                                                                     drucker_prager_parameters.max_yield_stress);

          // Step 4b: select if yield viscosity is based on Drucker Prager or stress limiter rheology
          double viscosity_yield = viscosity_pre_yield;
          switch (yield_type)
            {
              case stress_limiter:
              {
                //Step 4b-1: always rescale the viscosity back to the yield surface
                const double viscosity_limiter = yield_stress / (2.0 * ref_strain_rate)
                                                 * std::pow((edot_ii/ref_strain_rate),
                                                            1./exponents_stress_limiter[j] - 1.0);
                viscosity_yield = 1. / ( 1./viscosity_limiter + 1./viscosity_pre_yield);
                break;
              }
              case drucker_prager:
              {
                // Step 4b-2: if the current stress is greater than the yield stress,
                // rescale the viscosity back to yield surface
                if (current_stress >= yield_stress)
                  {
                    viscosity_yield = drucker_prager_plasticity.compute_viscosity(current_cohesion,
                                                                                  current_friction,
                                                                                  std::max(in.pressure[i], 0.0),
                                                                                  current_edot_ii,
                                                                                  drucker_prager_parameters.max_yield_stress);
                    composition_yielding[j] = true;
                  }
                break;
              }
              default:
              {
                AssertThrow(false, ExcNotImplemented());
                break;
              }
            }

          // Step 5: limit the viscosity with specified minimum and maximum bounds
          composition_viscosities[j] = std::min(std::max(viscosity_yield, min_visc), max_visc);
        }
      return std::make_pair (composition_viscosities, composition_yielding);
    }


    template <int dim>
    void
    ViscoPlastic<dim>::
    fill_plastic_outputs(const unsigned int i,
                         const std::vector<double> &volume_fractions,
                         const bool plastic_yielding,
                         const MaterialModel::MaterialModelInputs<dim> &in,
                         MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      PlasticAdditionalOutputs<dim> *plastic_out = out.template get_additional_output<PlasticAdditionalOutputs<dim> >();

      if (plastic_out != nullptr)
        {
          plastic_out->cohesions[i] = 0;
          plastic_out->friction_angles[i] = 0;
          plastic_out->yielding[i] = plastic_yielding ? 1 : 0;

          // set to weakened values, or unweakened values when strain weakening is not used
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            {
              // Calculate the strain weakening factors and weakened values
              const std::array<double, 3> weakening_factors = strain_rheology.compute_strain_weakening_factors(j, in.composition[i]);
              plastic_out->cohesions[i]   += volume_fractions[j] * (drucker_prager_parameters.cohesions[j] * weakening_factors[0]);
              // Also convert radians to degrees
              plastic_out->friction_angles[i] += 180.0/numbers::PI * volume_fractions[j] * (drucker_prager_parameters.angles_internal_friction[j] * weakening_factors[1]);
            }
        }
    }



    template <int dim>
    void
    ViscoPlastic<dim>::
    compute_viscosity_derivatives(const unsigned int i,
                                  const std::vector<double> &volume_fractions,
                                  const std::vector<double> &composition_viscosities,
                                  const MaterialModel::MaterialModelInputs<dim> &in,
                                  MaterialModel::MaterialModelOutputs<dim> &out,
                                  const std::vector<double> &phase_function_values) const
    {
      MaterialModel::MaterialModelDerivatives<dim> *derivatives =
        out.template get_additional_output<MaterialModel::MaterialModelDerivatives<dim> >();

      if (derivatives != nullptr)
        {
          // compute derivatives if necessary
          std::vector<SymmetricTensor<2,dim> > composition_viscosities_derivatives(volume_fractions.size());
          std::vector<double> composition_dviscosities_dpressure(volume_fractions.size());

          const double finite_difference_accuracy = 1e-7;

          // A new material model inputs variable that uses the strain rate and pressure difference.
          MaterialModel::MaterialModelInputs<dim> in_derivatives = in;

          // For each independent component, compute the derivative.
          for (unsigned int component = 0; component < SymmetricTensor<2,dim>::n_independent_components; ++component)
            {
              const TableIndices<2> strain_rate_indices = SymmetricTensor<2,dim>::unrolled_to_component_indices (component);

              // components that are not on the diagonal are multiplied by 0.5, because the symmetric tensor
              // is modified by 0.5 in both symmetric directions (xy/yx) simultaneously and we compute the combined
              // derivative
              const SymmetricTensor<2,dim> strain_rate_difference = in.strain_rate[i]
                                                                    + std::max(std::fabs(in.strain_rate[i][strain_rate_indices]), min_strain_rate)
                                                                    * (component > dim-1 ? 0.5 : 1 )
                                                                    * finite_difference_accuracy
                                                                    * Utilities::nth_basis_for_symmetric_tensors<dim>(component);

              in_derivatives.strain_rate[i] = strain_rate_difference;

              std::vector<double> eta_component =
                calculate_isostrain_viscosities(in_derivatives, i, volume_fractions,
                                                viscous_flow_law, yield_mechanism,
                                                phase_function_values).first;

              // For each composition of the independent component, compute the derivative.
              for (unsigned int composition_index = 0; composition_index < eta_component.size(); ++composition_index)
                {
                  // compute the difference between the viscosity with and without the strain-rate difference.
                  double viscosity_derivative = eta_component[composition_index] - composition_viscosities[composition_index];
                  if (viscosity_derivative != 0)
                    {
                      // when the difference is non-zero, divide by the difference.
                      viscosity_derivative /= std::max(std::fabs(strain_rate_difference[strain_rate_indices]), min_strain_rate)
                                              * finite_difference_accuracy;
                    }
                  composition_viscosities_derivatives[composition_index][strain_rate_indices] = viscosity_derivative;
                }
            }

          /**
           * Now compute the derivative of the viscosity to the pressure
           */
          const double pressure_difference = in.pressure[i] + (std::fabs(in.pressure[i]) * finite_difference_accuracy);

          in_derivatives.pressure[i] = pressure_difference;

          // Modify the in_derivatives object again to take the original strain rate.
          in_derivatives.strain_rate[i] = in.strain_rate[i];

          const std::vector<double> viscosity_difference =
            calculate_isostrain_viscosities(in_derivatives, i, volume_fractions,
                                            viscous_flow_law, yield_mechanism,
                                            phase_function_values).first;

          for (unsigned int composition_index = 0; composition_index < viscosity_difference.size(); ++composition_index)
            {
              double viscosity_derivative = viscosity_difference[composition_index] - composition_viscosities[composition_index];
              if (viscosity_difference[composition_index] != 0)
                {
                  if (in.pressure[i] != 0)
                    {
                      viscosity_derivative /= std::fabs(in.pressure[i]) * finite_difference_accuracy;
                    }
                  else
                    {
                      viscosity_derivative = 0;
                    }
                }
              composition_dviscosities_dpressure[composition_index] = viscosity_derivative;
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
                                                             composition_viscosities,
                                                             composition_viscosities_derivatives,
                                                             viscosity_averaging_p);
          derivatives->viscosity_derivative_wrt_pressure[i] =
            Utilities::derivative_of_weighted_p_norm_average(out.viscosities[i],
                                                             volume_fractions,
                                                             composition_viscosities,
                                                             composition_dviscosities_dpressure,
                                                             viscosity_averaging_p);
        }
    }


    template <int dim>
    ComponentMask
    ViscoPlastic<dim>::
    get_volumetric_composition_mask() const
    {
      // Store which components to exclude during the volume fraction computation.
      ComponentMask composition_mask = strain_rheology.get_strain_composition_mask();

      if (use_elasticity)
        {
          for (unsigned int i = 0; i < SymmetricTensor<2,dim>::n_independent_components ; ++i)
            composition_mask.set(i,false);
        }

      return composition_mask;
    }



    template <int dim>
    void
    ViscoPlastic<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // Store which components do not represent volumetric compositions (e.g. strain components).
      const ComponentMask volumetric_compositions = get_volumetric_composition_mask();

      EquationOfStateOutputs<dim> eos_outputs (this->n_compositional_fields()+1);
      EquationOfStateOutputs<dim> eos_outputs_all_phases (this->n_compositional_fields()+1+phase_function.n_phase_transitions());

      std::vector<double> average_elastic_shear_moduli (in.n_evaluation_points());

      // Store value of phase function for each phase and composition
      // While the number of phases is fixed, the value of the phase function is updated for every point
      std::vector<double> phase_function_values(phase_function.n_phase_transitions(), 0.0);

      // Loop through all requested points
      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          // First compute the equation of state variables and thermodynamic properties
          equation_of_state.evaluate(in, i, eos_outputs_all_phases);

          const double gravity_norm = this->get_gravity_model().gravity_vector(in.position[i]).norm();
          const double reference_density = (this->get_adiabatic_conditions().is_initialized())
                                           ?
                                           this->get_adiabatic_conditions().density(in.position[i])
                                           :
                                           eos_outputs_all_phases.densities[0];

          // The phase index is set to invalid_unsigned_int, because it is only used internally
          // in phase_average_equation_of_state_outputs to loop over all existing phases
          MaterialUtilities::PhaseFunctionInputs<dim> phase_inputs(in.temperature[i],
                                                                   in.pressure[i],
                                                                   this->get_geometry_model().depth(in.position[i]),
                                                                   gravity_norm*reference_density,
                                                                   numbers::invalid_unsigned_int);

          // Compute value of phase functions
          for (unsigned int j=0; j < phase_function.n_phase_transitions(); j++)
            {
              phase_inputs.phase_index = j;
              phase_function_values[j] = phase_function.compute_value(phase_inputs);
            }

          // Average by value of gamma function to get value of compositions
          phase_average_equation_of_state_outputs(eos_outputs_all_phases,
                                                  phase_function_values,
                                                  phase_function.n_phase_transitions_for_each_composition(),
                                                  eos_outputs);

          const std::vector<double> volume_fractions = MaterialUtilities::compute_volume_fractions(in.composition[i], volumetric_compositions);

          // not strictly correct if thermal expansivities are different, since we are interpreting
          // these compositions as volume fractions, but the error introduced should not be too bad.
          out.densities[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.densities, MaterialUtilities::arithmetic);
          out.thermal_expansion_coefficients[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.thermal_expansion_coefficients, MaterialUtilities::arithmetic);
          out.specific_heat[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.specific_heat_capacities, MaterialUtilities::arithmetic);

          if (define_conductivities == false)
            {
              double thermal_diffusivity = 0.0;

              for (unsigned int j=0; j < volume_fractions.size(); ++j)
                thermal_diffusivity += volume_fractions[j] * thermal_diffusivities[j];

              // Thermal conductivity at the given positions. If the temperature equation uses
              // the reference density profile formulation, use the reference density to
              // calculate thermal conductivity. Otherwise, use the real density. If the adiabatic
              // conditions are not yet initialized, the real density will still be used.
              if (this->get_parameters().formulation_temperature_equation ==
                  Parameters<dim>::Formulation::TemperatureEquation::reference_density_profile &&
                  this->get_adiabatic_conditions().is_initialized())
                out.thermal_conductivities[i] = thermal_diffusivity * out.specific_heat[i] *
                                                this->get_adiabatic_conditions().density(in.position[i]);
              else
                out.thermal_conductivities[i] = thermal_diffusivity * out.specific_heat[i] * out.densities[i];
            }
          else
            {
              // Use thermal conductivity values specified in the parameter file, if this
              // option was selected.
              out.thermal_conductivities[i] = MaterialUtilities::average_value (volume_fractions, thermal_conductivities, MaterialUtilities::arithmetic);
            }

          out.compressibilities[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.compressibilities, MaterialUtilities::arithmetic);
          out.entropy_derivative_pressure[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.entropy_derivative_pressure, MaterialUtilities::arithmetic);
          out.entropy_derivative_temperature[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.entropy_derivative_temperature, MaterialUtilities::arithmetic);

          // Compute the effective viscosity if requested and retrieve whether the material is plastically yielding
          bool plastic_yielding = false;
          if (in.requests_property(MaterialProperties::viscosity))
            {
              // Currently, the viscosities for each of the compositional fields are calculated assuming
              // isostrain amongst all compositions, allowing calculation of the viscosity ratio.
              // TODO: This is only consistent with viscosity averaging if the arithmetic averaging
              // scheme is chosen. It would be useful to have a function to calculate isostress viscosities.
              const std::pair<std::vector<double>, std::vector<bool> > calculate_viscosities =
                calculate_isostrain_viscosities(in, i, volume_fractions, viscous_flow_law,
                                                yield_mechanism, phase_function_values);

              // The isostrain condition implies that the viscosity averaging should be arithmetic (see above).
              // We have given the user freedom to apply alternative bounds, because in diffusion-dominated
              // creep (where n_diff=1) viscosities are stress and strain-rate independent, so the calculation
              // of compositional field viscosities is consistent with any averaging scheme.
              out.viscosities[i] = MaterialUtilities::average_value(volume_fractions, calculate_viscosities.first, viscosity_averaging);

              // Decide based on the maximum composition if material is yielding.
              // This avoids for example division by zero for harmonic averaging (as plastic_yielding
              // holds values that are either 0 or 1), but might not be consistent with the viscosity
              // averaging chosen.
              std::vector<double>::const_iterator max_composition = std::max_element(volume_fractions.begin(),volume_fractions.end());
              plastic_yielding = calculate_viscosities.second[std::distance(volume_fractions.begin(),max_composition)];

              // Compute viscosity derivatives if they are requested
              if (MaterialModel::MaterialModelDerivatives<dim> *derivatives =
                    out.template get_additional_output<MaterialModel::MaterialModelDerivatives<dim> >())
                compute_viscosity_derivatives(i, volume_fractions, calculate_viscosities.first, in, out, phase_function_values);
            }

          // Now compute changes in the compositional fields (i.e. the accumulated strain).
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;

          // Calculate changes in strain invariants and update the reaction terms
          strain_rheology.fill_reaction_outputs(in, i, min_strain_rate, plastic_yielding, out);

          // Fill plastic outputs if they exist.
          fill_plastic_outputs(i,volume_fractions,plastic_yielding,in,out);

          if (use_elasticity)
            {
              // Compute average elastic shear modulus
              average_elastic_shear_moduli[i] = MaterialUtilities::average_value(volume_fractions,
                                                                                 elastic_rheology.get_elastic_shear_moduli(),
                                                                                 viscosity_averaging);

              // Fill the material properties that are part of the elastic additional outputs
              if (ElasticAdditionalOutputs<dim> *elastic_out = out.template get_additional_output<ElasticAdditionalOutputs<dim> >())
                {
                  elastic_out->elastic_shear_moduli[i] = average_elastic_shear_moduli[i];
                }
            }
        }

      // If we use the full strain tensor, compute the change in the individual tensor components.
      strain_rheology.compute_finite_strain_reaction_terms(in, out);

      if (use_elasticity)
        {
          elastic_rheology.fill_elastic_force_outputs(in, average_elastic_shear_moduli, out);
          elastic_rheology.fill_reaction_outputs(in, average_elastic_shear_moduli, out);
        }
    }

    template <int dim>
    double
    ViscoPlastic<dim>::
    reference_viscosity () const
    {
      return ref_visc;
    }

    template <int dim>
    bool
    ViscoPlastic<dim>::
    is_compressible () const
    {
      return equation_of_state.is_compressible();
    }

    template <int dim>
    double ViscoPlastic<dim>::
    get_min_strain_rate () const
    {
      return min_strain_rate;
    }

    template <int dim>
    void
    ViscoPlastic<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Visco Plastic");
        {
          MaterialUtilities::PhaseFunction<dim>::declare_parameters(prm);

          EquationOfState::MulticomponentIncompressible<dim>::declare_parameters (prm);

          Rheology::StrainDependent<dim>::declare_parameters (prm);

          Rheology::Elasticity<dim>::declare_parameters (prm);

          // Reference and minimum/maximum values
          prm.declare_entry ("Minimum strain rate", "1.0e-20", Patterns::Double (0.),
                             "Stabilizes strain dependent viscosity. Units: \\si{\\per\\second}.");
          prm.declare_entry ("Reference strain rate","1.0e-15",Patterns::Double (0.),
                             "Reference strain rate for first time step. Units: \\si{\\per\\second}.");
          prm.declare_entry ("Minimum viscosity", "1e17", Patterns::Double (0.),
                             "Lower cutoff for effective viscosity. Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Maximum viscosity", "1e28", Patterns::Double (0.),
                             "Upper cutoff for effective viscosity. Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Reference viscosity", "1e22", Patterns::Double (0.),
                             "Reference viscosity for nondimensionalization. "
                             "To understand how pressure scaling works, take a look at "
                             "\\cite{KHB12}. In particular, the value of this parameter "
                             "would not affect the solution computed by \\aspect{} if "
                             "we could do arithmetic exactly; however, computers do "
                             "arithmetic in finite precision, and consequently we need to "
                             "scale quantities in ways so that their magnitudes are "
                             "roughly the same. As explained in \\cite{KHB12}, we scale "
                             "the pressure during some computations (never visible by "
                             "users) by a factor that involves a reference viscosity. This "
                             "parameter describes this reference viscosity."
                             "\n\n"
                             "For problems with a constant viscosity, you will generally want "
                             "to choose the reference viscosity equal to the actual viscosity. "
                             "For problems with a variable viscosity, the reference viscosity "
                             "should be a value that adequately represents the order of "
                             "magnitude of the viscosities that appear, such as an average "
                             "value or the value one would use to compute a Rayleigh number."
                             "\n\n"
                             "Units: \\si{\\pascal\\second}.");

          // Equation of state parameters
          prm.declare_entry ("Thermal diffusivities", "0.8e-6",
                             Patterns::List(Patterns::Double (0.)),
                             "List of thermal diffusivities, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value.  "
                             "Units: \\si{\\meter\\squared\\per\\second}.");
          prm.declare_entry ("Define thermal conductivities","false",
                             Patterns::Bool (),
                             "Whether to directly define thermal conductivities for each compositional field "
                             "instead of calculating the values through the specified thermal diffusivities, "
                             "densities, and heat capacities. ");
          prm.declare_entry ("Thermal conductivities", "3.0",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal conductivities, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. "
                             "Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");

          // Rheological parameters
          prm.declare_entry ("Viscosity averaging scheme", "harmonic",
                             Patterns::Selection("arithmetic|harmonic|geometric|maximum composition"),
                             "When more than one compositional field is present at a point "
                             "with different viscosities, we need to come up with an average "
                             "viscosity at that point.  Select a weighted harmonic, arithmetic, "
                             "geometric, or maximum composition.");
          prm.declare_entry ("Viscous flow law", "composite",
                             Patterns::Selection("diffusion|dislocation|frank kamenetskii|composite"),
                             "Select what type of viscosity law to use between diffusion, "
                             "dislocation, frank kamenetskii, and composite options. Soon there will be an option "
                             "to select a specific flow law for each assigned composition ");
          prm.declare_entry ("Yield mechanism", "drucker",
                             Patterns::Selection("drucker|limiter"),
                             "Select what type of yield mechanism to use between Drucker Prager "
                             "and stress limiter options.");

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

          // Drucker Prager plasticity parameters
          Rheology::DruckerPrager<dim>::declare_parameters(prm);

          // Stress limiter parameters
          prm.declare_entry ("Stress limiter exponents", "1.0",
                             Patterns::List(Patterns::Double (0.)),
                             "List of stress limiter exponents, $n_{\\text{lim}}$, "
                             "for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "Units: none.");

          // Temperature in viscosity laws to include an adiabat (note units of K/Pa)
          prm.declare_entry ("Adiabat temperature gradient for viscosity", "0.0", Patterns::Double (0.),
                             "Add an adiabatic temperature gradient to the temperature used in the flow law "
                             "so that the activation volume is consistent with what one would use in a "
                             "earth-like (compressible) model. Default is set so this is off. "
                             "Note that this is a linear approximation of the real adiabatic gradient, which "
                             "is okay for the upper mantle, but is not really accurate for the lower mantle. "
                             "Using a pressure gradient of 32436 Pa/m, then a value of "
                             "0.3 K/km = 0.0003 K/m = 9.24e-09 K/Pa gives an earth-like adiabat."
                             "Units: \\si{\\kelvin\\per\\pascal}.");

          prm.declare_entry ("Include viscoelasticity", "false",
                             Patterns::Bool (),
                             "Whether to include elastic effects in the rheological formulation.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    ViscoPlastic<dim>::parse_parameters (ParameterHandler &prm)
    {
      // increment by one for background:
      const unsigned int n_fields = this->n_compositional_fields() + 1;

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Visco Plastic");
        {
          // Phase transition parameters
          phase_function.initialize_simulator (this->get_simulator());
          phase_function.parse_parameters (prm);

          std::vector<unsigned int> n_phase_transitions_for_each_composition
          (phase_function.n_phase_transitions_for_each_composition());

          // We require one more entry for density, etc as there are phase transitions
          // (for the low-pressure phase before any transition).
          for (unsigned int i=0; i<n_phase_transitions_for_each_composition.size(); ++i)
            n_phase_transitions_for_each_composition[i] += 1;

          // Equation of state parameters
          equation_of_state.initialize_simulator (this->get_simulator());
          equation_of_state.parse_parameters (prm,
                                              std::make_shared<std::vector<unsigned int>>(n_phase_transitions_for_each_composition));

          strain_rheology.initialize_simulator (this->get_simulator());
          strain_rheology.parse_parameters(prm);

          use_elasticity = prm.get_bool ("Include viscoelasticity");

          if (use_elasticity)
            {
              elastic_rheology.initialize_simulator (this->get_simulator());
              elastic_rheology.parse_parameters(prm);
            }

          // Reference and minimum/maximum values
          min_strain_rate = prm.get_double("Minimum strain rate");
          ref_strain_rate = prm.get_double("Reference strain rate");
          min_visc = prm.get_double ("Minimum viscosity");
          max_visc = prm.get_double ("Maximum viscosity");
          ref_visc = prm.get_double ("Reference viscosity");

          thermal_diffusivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal diffusivities"))),
                                                                          n_fields,
                                                                          "Thermal diffusivities");

          define_conductivities = prm.get_bool ("Define thermal conductivities");

          thermal_conductivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal conductivities"))),
                                                                           n_fields,
                                                                           "Thermal conductivities");


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
          else
            AssertThrow(false, ExcMessage("Not a valid viscous flow law"));

          // Rheological parameters
          if (prm.get ("Yield mechanism") == "drucker")
            yield_mechanism = drucker_prager;
          else if (prm.get ("Yield mechanism") == "limiter")
            yield_mechanism = stress_limiter;
          else
            AssertThrow(false, ExcMessage("Not a valid yield mechanism."));

          AssertThrow(use_elasticity == false || yield_mechanism == drucker_prager,
                      ExcMessage("Elastic behavior is only tested with the "
                                 "'drucker prager' plasticity option."));

          // Rheological parameters
          // Diffusion creep parameters
          diffusion_creep.initialize_simulator (this->get_simulator());
          diffusion_creep.parse_parameters(prm, std::make_shared<std::vector<unsigned int>>(n_phase_transitions_for_each_composition));

          // Dislocation creep parameters
          dislocation_creep.initialize_simulator (this->get_simulator());
          dislocation_creep.parse_parameters(prm, std::make_shared<std::vector<unsigned int>>(n_phase_transitions_for_each_composition));

          // Frank Kamenetskii viscosity parameters
          if (viscous_flow_law == frank_kamenetskii)
            {
              frank_kamenetskii_rheology = std_cxx14::make_unique<Rheology::FrankKamenetskii<dim>>();
              frank_kamenetskii_rheology->initialize_simulator (this->get_simulator());
              frank_kamenetskii_rheology->parse_parameters(prm);
            }

          // Peierls creep parameters
          use_peierls_creep = prm.get_bool ("Include Peierls creep");
          if (use_peierls_creep)
            {
              peierls_creep = std_cxx14::make_unique<Rheology::PeierlsCreep<dim>>();
              peierls_creep->initialize_simulator (this->get_simulator());
              peierls_creep->parse_parameters(prm);
            }

          // Constant viscosity prefactor parameters
          constant_viscosity_prefactors.initialize_simulator (this->get_simulator());
          constant_viscosity_prefactors.parse_parameters(prm);

          // Plasticity parameters
          drucker_prager_parameters = drucker_prager_plasticity.parse_parameters(this->n_compositional_fields()+1,
                                                                                 prm);

          // Stress limiter parameter
          exponents_stress_limiter  = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Stress limiter exponents"))),
                                                                              n_fields,
                                                                              "Stress limiter exponents");

          // Include an adiabat temperature gradient in flow laws
          adiabatic_temperature_gradient_for_viscosity = prm.get_double("Adiabat temperature gradient for viscosity");
          if (this->get_heating_model_manager().adiabatic_heating_enabled())
            AssertThrow (adiabatic_temperature_gradient_for_viscosity == 0.0,
                         ExcMessage("If adiabatic heating is enabled you should not add another adiabatic gradient"
                                    "to the temperature for computing the viscosity, because the ambient"
                                    "temperature profile already includes the adiabatic gradient."));


        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::strain_rate | NonlinearDependence::compositional_fields;
      this->model_dependence.density = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
    }

    template <int dim>
    void
    ViscoPlastic<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output<PlasticAdditionalOutputs<dim> >() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::PlasticAdditionalOutputs<dim>> (n_points));
        }
      if (use_elasticity)
        elastic_rheology.create_elastic_outputs(out);
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(ViscoPlastic,
                                   "visco plastic",
                                   "An implementation of an incompressible visco(elastic)-plastic rheology "
                                   "with options for selecting dislocation creep, diffusion creep or "
                                   "composite viscous flow laws. Prior to yielding, one may select to "
                                   "modify the viscosity to account for viscoelastic effects. Plasticity "
                                   "limits viscous stresses through a Drucker Prager yield criterion. "
                                   "Note that this material model is based heavily on the "
                                   "DiffusionDislocation (Bob Myhill), DruckerPrager "
                                   "(Anne Glerum), and Viscoelastic (John Naliboff) material models. "
                                   "\n\n "
                                   "The viscosity for dislocation or diffusion creep is defined as "
                                   "\\[v = \\frac 12 A^{-\\frac{1}{n}} d^{\\frac{m}{n}} "
                                   "\\dot{\\varepsilon}_{ii}^{\\frac{1-n}{n}} "
                                   "\\exp\\left(\\frac{E + PV}{nRT}\\right)\\] "
                                   "where $A$ is the prefactor, $n$ is the stress exponent, "
                                   "$\\dot{\\varepsilon}_{ii}$ is the square root of the deviatoric "
                                   "strain rate tensor second invariant, $d$ is grain size, "
                                   "$m$ is the grain size exponent, $E$ is activation energy, "
                                   "$V$ is activation volume, $P$ is pressure, $R$ is the gas "
                                   "exponent and $T$ is temperature. "
                                   "This form of the viscosity equation is commonly used in "
                                   "geodynamic simulations. See, for example, Billen and Hirth "
                                   "(2007), G3, 8, Q08012. Significantly, other studies may use "
                                   "slightly different forms of the viscosity equation leading to "
                                   "variations in how specific terms are defined or combined. For "
                                   "example, the grain size exponent should always be positive in "
                                   "the diffusion viscosity equation used here, while other studies "
                                   "place the grain size term in the denominator and invert the sign "
                                   "of the grain size exponent. When examining previous work, one "
                                   "should carefully check how the viscous prefactor and grain size "
                                   "terms are defined. "
                                   "\n\n "
                                   "One may select to use the diffusion ($v_{\\text{diff}}$; $n=1$, $m!=0$), "
                                   "dislocation ($v_{\\text{disl}}$, $n>1$, $m=0$) or composite "
                                   "$\\frac{v_{\\text{diff}} v_{\\text{disl}}}{v_{\\text{diff}}+v_{\\text{disl}}}$ equation form. "
                                   "\n\n "
                                   "The diffusion and dislocation prefactors can be weakened with a factor "
                                   "between 0 and 1 according to the total or the viscous strain only. "
                                   "\n\n "
                                   "Viscosity is limited through one of two different `yielding' mechanisms. "
                                   "\n\n"
                                   "The first plasticity mechanism limits viscous stress through a "
                                   "Drucker Prager yield criterion, where the yield stress in 3D is  "
                                   "$\\sigma_y = \\frac{6C\\cos(\\phi) + 2P\\sin(\\phi)} "
                                   "{\\sqrt(3)(3+\\sin(\\phi))}$ "
                                   "and "
                                   "$\\sigma_y = C\\cos(\\phi) + P\\sin(\\phi)$ "
                                   "in 2D. Above, $C$ is cohesion and $\\phi$  is the angle of "
                                   "internal friction.  Note that the 2D form is equivalent to the "
                                   "Mohr Coulomb yield surface.  If $\\phi$ is 0, the yield stress "
                                   "is fixed and equal to the cohesion (Von Mises yield criterion). "
                                   "When the viscous stress ($2v{\\varepsilon}_{ii}$) exceeds "
                                   "the yield stress, the viscosity is rescaled back to the yield "
                                   "surface: $v_{y}=\\sigma_{y}/(2{\\varepsilon}_{ii})$. "
                                   "This form of plasticity is commonly used in geodynamic models. "
                                   "See, for example, Thieulot, C. (2011), PEPI 188, pp. 47-68. "
                                   "\n\n"
                                   "The user has the option to linearly reduce the cohesion and "
                                   "internal friction angle as a function of the finite strain magnitude. "
                                   "The finite strain invariant or full strain tensor is calculated through "
                                   "compositional fields within the material model. This implementation is "
                                   "identical to the compositional field finite strain plugin and cookbook "
                                   "described in the manual (author: Gassmoeller, Dannberg). If the user selects to track "
                                   "the finite strain invariant ($e_{ii}$), a single compositional field tracks "
                                   "the value derived from $e_{ii}^t = (e_{ii})^{(t-1)} + \\dot{e}_{ii}\\; dt$, where $t$ and $t-1$ "
                                   "are the current and prior time steps, $\\dot{e}_{ii}$ is the second invariant of the "
                                   "strain rate tensor and $dt$ is the time step size. In the case of the "
                                   "full strain tensor $F$, the finite strain magnitude is derived from the "
                                   "second invariant of the symmetric stretching tensor $L$, where "
                                   "$L = F [F]^T$. The user must specify a single compositional "
                                   "field for the finite strain invariant or multiple fields (4 in 2D, 9 in 3D) "
                                   "for the finite strain tensor. These field(s) must be the first listed "
                                   "compositional fields in the parameter file. Note that one or more of the finite strain "
                                   "tensor components must be assigned a non-zero value initially. This value can be "
                                   "be quite small (e.g., 1.e-8), but still non-zero. While the option to track and use "
                                   "the full finite strain tensor exists, tracking the associated compositional fields "
                                   "is computationally expensive in 3D. Similarly, the finite strain magnitudes "
                                   "may in fact decrease if the orientation of the deformation field switches "
                                   "through time. Consequently, the ideal solution is track the finite strain "
                                   "invariant (single compositional) field within the material and track "
                                   "the full finite strain tensor through particles."
                                   "When only the second invariant of the strain is tracked, one has the option to "
                                   "track the full strain or only the plastic strain. In the latter case, strain is only tracked "
                                   "in case the material is plastically yielding, i.e. the viscous stress > yield stress. "
                                   "\n\n"
                                   "Viscous stress may also be limited by a non-linear stress limiter "
                                   "that has a form similar to the Peierls creep mechanism. "
                                   "This stress limiter assigns an effective viscosity "
                                   "$\\sigma_{\\text{eff}} = \\frac{\\tau_y}{2\\varepsilon_y} "
                                   "{\\frac{\\varepsilon_{ii}}{\\varepsilon_y}}^{\\frac{1}{n_y}-1}$ "
                                   "Above $\\tau_y$ is a yield stress, $\\varepsilon_y$ is the "
                                   "reference strain rate, $\\varepsilon_{ii}$ is the strain rate "
                                   "and $n_y$ is the stress limiter exponent.  The yield stress, "
                                   "$\\tau_y$, is defined through the Drucker Prager yield criterion "
                                   "formulation. This method of limiting viscous stress has been used "
                                   "in various forms within the geodynamic literature \\cite{chri92,vavv02,cibi13,cibi15}."
                                   "When $n_y$ is 1, it essentially becomes a linear viscosity model, "
                                   "and in the limit $n_y\\rightarrow \\infty$ it converges to the "
                                   "standard viscosity rescaling method (concretely, values $n_y>20$ "
                                   "are large enough)."
                                   "\n\n "
                                   "The visco-plastic rheology described above may also be modified to include "
                                   "viscoelastic deformation, thus producing a viscoelastic plastic constitutive "
                                   "relationship. "
                                   "\n\n "
                                   "The viscoelastic rheology behavior takes into account the elastic shear "
                                   "strength (e.g., shear modulus), while the tensile and volumetric "
                                   "strength (e.g., Young's and bulk modulus) are not considered. The "
                                   "model is incompressible and allows specifying an arbitrary number "
                                   "of compositional fields, where each field represents a different "
                                   "rock type or component of the viscoelastic stress tensor. The stress "
                                   "tensor in 2D and 3D, respectively, contains 3 or 6 components. The "
                                   "compositional fields representing these components must be named "
                                   "and listed in a very specific format, which is designed to minimize "
                                   "mislabeling stress tensor components as distinct 'compositional "
                                   "rock types' (or vice versa). For 2D models, the first three "
                                   "compositional fields must be labeled 'stress\\_xx', 'stress\\_yy' and 'stress\\_xy'. "
                                   "In 3D, the first six compositional fields must be labeled 'stress\\_xx', "
                                   "'stress\\_yy', 'stress\\_zz', 'stress\\_xy', 'stress\\_xz', 'stress\\_yz'. "
                                   "\n\n "
                                   "Combining this viscoelasticity implementation with non-linear viscous flow "
                                   "and plasticity produces a constitutive relationship commonly referred to "
                                   "as partial elastoviscoplastic (e.g., pEVP) in the geodynamics community. "
                                   "While extensively discussed and applied within the geodynamics "
                                   "literature, notable references include: "
                                   "Moresi et al. (2003), J. Comp. Phys., v. 184, p. 476-497. "
                                   "Gerya and Yuen (2007), Phys. Earth. Planet. Inter., v. 163, p. 83-105. "
                                   "Gerya (2010), Introduction to Numerical Geodynamic Modeling. "
                                   "Kaus (2010), Tectonophysics, v. 484, p. 36-47. "
                                   "Choi et al. (2013), J. Geophys. Res., v. 118, p. 2429-2444. "
                                   "Keller et al. (2013), Geophys. J. Int., v. 195, p. 1406-1442. "
                                   "\n\n "
                                   "The overview below directly follows Moresi et al. (2003) eqns. 23-38. "
                                   "However, an important distinction between this material model and "
                                   "the studies above is the use of compositional fields, rather than "
                                   "tracers, to track individual components of the viscoelastic stress "
                                   "tensor. The material model will be updated when an option to track "
                                   "and calculate viscoelastic stresses with tracers is implemented. "
                                   "\n\n "
                                   "Moresi et al. (2003) begins (eqn. 23) by writing the deviatoric "
                                   "rate of deformation ($\\hat{D}$) as the sum of elastic "
                                   "($\\hat{D_{e}}$) and viscous ($\\hat{D_{v}}$) components: "
                                   "$\\hat{D} = \\hat{D_{e}} + \\hat{D_{v}}$.  "
                                   "These terms further decompose into "
                                   "$\\hat{D_{v}} = \\frac{\\tau}{2\\eta}$ and "
                                   "$\\hat{D_{e}} = \\frac{\\overset{\\nabla}{\\tau}}{2\\mu}$, where "
                                   "$\\tau$ is the viscous deviatoric stress, $\\eta$ is the shear viscosity, "
                                   "$\\mu$ is the shear modulus and $\\overset{\\nabla}{\\tau}$ is the "
                                   "Jaumann corotational stress rate. This later term (eqn. 24) contains the "
                                   "time derivative of the deviatoric stress ($\\dot{\\tau}$) and terms that "
                                   "account for material spin (e.g., rotation) due to advection: "
                                   "$\\overset{\\nabla}{\\tau} = \\dot{\\tau} + {\\tau}W -W\\tau$. "
                                   "Above, $W$ is the material spin tensor (eqn. 25): "
                                   "$W_{ij} = \\frac{1}{2} \\left (\\frac{\\partial V_{i}}{\\partial x_{j}} - "
                                   "\\frac{\\partial V_{j}}{\\partial x_{i}} \\right )$. "
                                   "\n\n "
                                   "If plasticity is included, the deviatoric rate of deformation may be written as: "
                                   "$\\hat{D} = \\hat{D_{e}} + \\hat{D_{v}} + \\hat{D_{p}}$, where $\\hat{D_{p}}$ "
                                   "is the plastic component. $\\hat{D_{p}}$ decomposes to $\\frac{\\tau_{y}}{2\\eta_{y}}$, "
                                   "where $\\tau_{y}$ is the yield stress and $\\eta_{y}$ is the viscosity rescaled "
                                   "to the yield surface. "
                                   "The Jaumann stress-rate can also be approximated using terms from the "
                                   "previous time step ($t$) and current time step ($t + \\Delta t^{e}$): "
                                   "$\\smash[t]{\\overset{\\nabla}{\\tau}}^{t + \\Delta t^{e}} \\approx "
                                   "\\frac{\\tau^{t + \\Delta t^{e} - \\tau^{t}}}{\\Delta t^{e}} - "
                                   "W^{t}\\tau^{t} + \\tau^{t}W^{t}$. "
                                   "In this material model, the size of the time step above ($\\Delta t^{e}$) "
                                   "can be specified as the numerical time step size or an independent fixed time "
                                   "step. If the latter case is selected, the user has an option to apply a "
                                   "stress averaging scheme to account for the differences between the numerical "
                                   "and fixed elastic time step (eqn. 32). If one selects to use a fixed elastic time "
                                   "step throughout the model run, this can still be achieved by using CFL and "
                                   "maximum time step values that restrict the numerical time step to a specific time."
                                   "\n\n "
                                   "The formulation above allows rewriting the total rate of deformation (eqn. 29) as\n "
                                   "$\\tau^{t + \\Delta t^{e}} = \\eta_{eff} \\left ( "
                                   "2\\hat{D}^{t + \\triangle t^{e}} + \\frac{\\tau^{t}}{\\mu \\Delta t^{e}} + "
                                   "\\frac{W^{t}\\tau^{t} - \\tau^{t}W^{t}}{\\mu}  \\right )$. "
                                   "\n\n "
                                   "The effective viscosity (eqn. 28) is a function of the viscosity ($\\eta$), "
                                   "elastic time step size ($\\Delta t^{e}$) and shear relaxation time "
                                   "($ \\alpha = \\frac{\\eta}{\\mu} $): "
                                   "$\\eta_{eff} = \\eta \\frac{\\Delta t^{e}}{\\Delta t^{e} + \\alpha}$ "
                                   "The magnitude of the shear modulus thus controls how much the effective "
                                   "viscosity is reduced relative to the initial viscosity. "
                                   "\n\n "
                                   "Elastic effects are introduced into the governing Stokes equations through "
                                   "an elastic force term (eqn. 30) using stresses from the previous time step: "
                                   "$F^{e,t} = -\\frac{\\eta_{eff}}{\\mu \\Delta t^{e}} \\tau^{t}$. "
                                   "This force term is added onto the right-hand side force vector in the "
                                   "system of equations. "
                                   "\n\n "
                                   "When plastic yielding occurs, the effective viscosity in equation 29 and 30 is the "
                                   "plastic viscosity (equation 36). If the current stress is below the plastic "
                                   "yield stress, the effective viscosity is still as defined in equation 28. "
                                   "During non-linear iterations, we define the current stress prior to yielding "
                                   "(e.g., value compared to yield stress) as "
                                   "$\\tau^{t + \\Delta t^{e}} = \\eta_{eff} \\left ( 2\\hat{D}^{t + \\triangle t^{e}} + "
                                   "\\frac{\\tau^{t}}{\\mu \\Delta t^{e}} \\right ) $"
                                   "\n\n "
                                   "Compositional fields can each be assigned individual values of "
                                   "thermal diffusivity, heat capacity, density, thermal "
                                   "expansivity and rheological parameters. "
                                   "\n\n "
                                   "If more than one compositional field is present at a given "
                                   "point, viscosities are averaged with an arithmetic, geometric "
                                   "harmonic (default) or maximum composition scheme. "
                                   "\n\n "
                                   "The value for the components of this formula and additional "
                                   "parameters are read from the parameter file in subsection "
                                   " 'Material model/Visco Plastic'.")
  }
}
