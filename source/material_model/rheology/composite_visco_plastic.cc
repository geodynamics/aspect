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


#include <aspect/material_model/rheology/composite_visco_plastic.h>
#include <aspect/heating_model/shear_heating.h>
#include <aspect/material_model/utilities.h>
#include <aspect/utilities.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      // The following function is not yet used by ASPECT, so
      // it is commented out to avoid tester issues.
      //
      // It will provide useful output when a material model which
      // uses the CompositeViscoPlastic rheology has been implemented.
      //
      // The composite visco plastic rheology calculates the decomposed strain
      // rates for each of the following deformation mechanisms:
      // diffusion creep, dislocation creep, Peierls creep,
      // Drucker-Prager plasticity, Kelvin (damped) elasticity
      // and a constant (high) viscosity limiter.
      // The values are provided in this order as a vector of additional
      // outputs. If the user declares one or more mechanisms inactive
      // (by assigning use_mechanism = False) then the corresponding
      // strain rate output will be equal to zero.
      //namespace
      //{
      //  std::vector<std::string> make_strain_rate_additional_outputs_names()
      //  {
      //    std::vector<std::string> names;
      //    names.emplace_back("edot_diffusion");
      //    names.emplace_back("edot_dislocation");
      //    names.emplace_back("edot_peierls");
      //    names.emplace_back("edot_drucker_prager");
      //    names.emplace_back("edot_kelvin");
      //    names.emplace_back("edot_limiter");
      //    return names;
      //  }
      //}



      namespace ViscosityAveraging
      {
        ViscosityAveraging::Kind
        parse (const std::string &parameter_name,
               const ParameterHandler &prm)
        {
          ViscosityAveraging::Kind scheme;
          if (prm.get (parameter_name) == "isostress")
            scheme = isostress;
          else if (prm.get (parameter_name) == "isostrain")
            scheme = isostrain;
          else
            {
              AssertThrow(false, ExcMessage("Not a valid viscosity averaging scheme. Choose either isostress or isostrain."));
            }

          return scheme;
        }
      }



      template <int dim>
      CompositeViscoPlastic<dim>::CompositeViscoPlastic ()
        = default;



      template <int dim>
      double
      CompositeViscoPlastic<dim>::compute_inverse_kelvin_viscosity(const std::vector<double> &volume_fractions) const
      {
        double inverse_kelvin_viscosity = 0.;
        if (this->get_parameters().enable_elasticity)
          {
            // Take the volume-weighted harmonic average of the individual component
            // shear moduli, as required for isostress (Reuss) material averaging.
            const std::vector<double> &elastic_shear_moduli = elasticity->get_elastic_shear_moduli();
            const double elastic_shear_modulus = MaterialUtilities::average_value(volume_fractions, elastic_shear_moduli, MaterialUtilities::harmonic);
            inverse_kelvin_viscosity = 1./elasticity->calculate_elastic_viscosity(elastic_shear_modulus);
          }
        return inverse_kelvin_viscosity;
      }



      template <int dim>
      SymmetricTensor<2,dim>
      CompositeViscoPlastic<dim>::compute_effective_strain_rate(const SymmetricTensor<2,dim> &strain_rate,
                                                                const SymmetricTensor<2,dim> &elastic_stress,
                                                                const double inverse_kelvin_viscosity) const
      {
        return strain_rate + (0.5 * elastic_stress * inverse_kelvin_viscosity);
      }



      template <int dim>
      double
      CompositeViscoPlastic<dim>::compute_viscosity (const double pressure,
                                                     const double temperature,
                                                     const double grain_size,
                                                     const std::vector<double> &volume_fractions,
                                                     const SymmetricTensor<2,dim> &effective_strain_rate,
                                                     const double inverse_kelvin_viscosity,
                                                     std::vector<double> &partial_strain_rates,
                                                     const std::vector<double> &phase_function_values,
                                                     const std::vector<unsigned int> &n_phase_transitions_per_composition) const
      {
        double viscosity = 0.;

        // Make sure partial_strain_rates is filled with zeros and is the right length
        partial_strain_rates.resize(n_decomposed_strain_rates);
        std::fill(partial_strain_rates.begin(), partial_strain_rates.end(), 0);

        // Compute the viscosity and the partial strain rates
        // according to the isostress or isostrain viscosity averaging scheme.
        switch (viscosity_averaging_scheme)
          {
            case ViscosityAveraging::Kind::isostress:
            {
              viscosity = compute_isostress_viscosity (pressure,
                                                       temperature,
                                                       grain_size,
                                                       volume_fractions,
                                                       effective_strain_rate,
                                                       inverse_kelvin_viscosity,
                                                       partial_strain_rates,
                                                       phase_function_values,
                                                       n_phase_transitions_per_composition);
              break;
            }
            case ViscosityAveraging::Kind::isostrain:
            {
              viscosity = compute_isostrain_viscosity (pressure,
                                                       temperature,
                                                       grain_size,
                                                       volume_fractions,
                                                       effective_strain_rate,
                                                       partial_strain_rates,
                                                       phase_function_values,
                                                       n_phase_transitions_per_composition);
              break;
            }
            default:
            {
              Assert (false, ExcNotImplemented());
            }
          }
        return viscosity;
      }



      template <int dim>
      double
      CompositeViscoPlastic<dim>::compute_isostress_viscosity (const double pressure,
                                                               const double temperature,
                                                               const double grain_size,
                                                               const std::vector<double> &volume_fractions,
                                                               const SymmetricTensor<2,dim> &effective_strain_rate,
                                                               const double inverse_kelvin_viscosity,
                                                               std::vector<double> &partial_strain_rates,
                                                               const std::vector<double> &phase_function_values,
                                                               const std::vector<unsigned int> &n_phase_transitions_per_composition) const
      {
        // If strain rate is zero (like during the first time step) set it to some very small number
        // to prevent a division-by-zero, and a floating point exception.
        // Otherwise, calculate the square-root of the norm of the second invariant of the deviatoric-
        // strain rate (often simplified as epsilondot_ii)
        const double edot_ii = std::max(std::sqrt(std::max(-second_invariant(deviator(effective_strain_rate)), 0.)),
                                        minimum_strain_rate);
        const double log_edot_ii = std::log(edot_ii);

        std::vector<Rheology::DiffusionCreepParameters> diffusion_creep_parameters;
        std::vector<Rheology::DislocationCreepParameters> dislocation_creep_parameters;
        std::vector<Rheology::PeierlsCreepParameters> peierls_creep_parameters;
        std::vector<Rheology::DruckerPragerParameters> drucker_prager_parameters;

        // 1) Estimate the stress running through the creep elements and the
        // maximum viscosity element, whose viscosity is defined in parse_parameters.
        // The creep elements for each composition are arranged in series.
        // Taking the minimum viscosity from
        // all of these elements provides an excellent first approximation
        // to the true viscosity of that composition.
        // For a starting guess, we further assume that each composition
        // experiences the same strain rate.

        // The stress can then be calculated as 2 * viscoplastic_viscosity_guess * edot_ii
        double inverse_viscoplastic_viscosity_guess = 1./maximum_viscosity;

        double total_volume_fraction = 1.;
        std::vector<unsigned int> active_compositions;
        for (unsigned int composition = 0; composition < number_of_chemical_compositions; ++composition)
          {
            // Only include the contribution to the viscosity
            // from a given composition if the volume fraction exceeds
            // a certain (small) fraction.
            double viscoplastic_viscosity_guess_c = maximum_viscosity;
            const double volume_fraction = volume_fractions[composition];
            if (volume_fraction > 2. * std::numeric_limits<double>::epsilon())
              {
                active_compositions.push_back(composition);

                if (use_diffusion_creep)
                  {
                    diffusion_creep_parameters.push_back(diffusion_creep->compute_creep_parameters(composition, phase_function_values, n_phase_transitions_per_composition));
                    viscoplastic_viscosity_guess_c = std::min(diffusion_creep->compute_viscosity(pressure, temperature, grain_size, composition, phase_function_values, n_phase_transitions_per_composition), viscoplastic_viscosity_guess_c);
                  }

                if (use_dislocation_creep)
                  {
                    dislocation_creep_parameters.push_back(dislocation_creep->compute_creep_parameters(composition, phase_function_values, n_phase_transitions_per_composition));
                    viscoplastic_viscosity_guess_c = std::min(dislocation_creep->compute_viscosity(edot_ii, pressure, temperature, composition, phase_function_values, n_phase_transitions_per_composition), viscoplastic_viscosity_guess_c);
                  }

                if (use_peierls_creep)
                  {
                    peierls_creep_parameters.push_back(peierls_creep->compute_creep_parameters(composition));
                    viscoplastic_viscosity_guess_c = std::min(peierls_creep->compute_approximate_viscosity(edot_ii, pressure, temperature, composition), viscoplastic_viscosity_guess_c);
                  }

                if (use_drucker_prager)
                  {
                    Rheology::DruckerPragerParameters drucker_prager_parameters_c = drucker_prager->compute_drucker_prager_parameters(composition, phase_function_values, n_phase_transitions_per_composition);
                    drucker_prager_parameters.push_back(drucker_prager_parameters_c);
                    viscoplastic_viscosity_guess_c = std::min(drucker_prager->compute_viscosity(drucker_prager_parameters_c.cohesion,
                                                                                                drucker_prager_parameters_c.angle_internal_friction, pressure, edot_ii, drucker_prager_parameters_c.max_yield_stress),
                                                              viscoplastic_viscosity_guess_c);
                  }
                inverse_viscoplastic_viscosity_guess += volume_fraction / viscoplastic_viscosity_guess_c;
              }
            else
              {
                total_volume_fraction -= volume_fractions[composition];
              }
          }

        const double viscoplastic_viscosity_guess = 1. / inverse_viscoplastic_viscosity_guess;
        double viscoplastic_stress = 2. * viscoplastic_viscosity_guess * edot_ii / total_volume_fraction;
        double log_viscoplastic_stress = std::log(viscoplastic_stress);

        // 2) Calculate the log strain rate and first derivative with respect to
        // log stress for each element using the guessed creep stress.
        std::vector<std::array<std::pair<double, double>, 4>> log_edot_and_deriv(active_compositions.size());

        for (unsigned int i = 0; i < active_compositions.size(); ++i)
          {
            // In the isostress viscosity averaging formalism,
            // the total strain rate is the sum of the strain rates for each
            // composition i multiplied by the volume fraction of that
            // composition. The strain rate for each composition is the sum
            // of strain rates for each deformation mechanism j:
            // edot_total = sum_i (volume_fraction_i * (sum_j (edot_ij)))
            // This can be refactored:
            // edot_total = sum_i (sum_j (volume_fraction_i * (edot_ij)))
            // Therefore the logarithm of the contribution of each
            // component-deformation mechanism to the total strain rate
            // is:
            // ln(edot_total_contrib_ij) = ln(volume_fraction_i) + ln(edot_ij)

            const double log_volume_fraction = std::log(volume_fractions[active_compositions[i]]);

            if (use_diffusion_creep)
              {
                log_edot_and_deriv[i][0] = diffusion_creep->compute_log_strain_rate_and_derivative(log_viscoplastic_stress, pressure, temperature, grain_size, diffusion_creep_parameters[i]);
                log_edot_and_deriv[i][0].first += log_volume_fraction;
              }
            if (use_dislocation_creep)
              {
                log_edot_and_deriv[i][1] = dislocation_creep->compute_log_strain_rate_and_derivative(log_viscoplastic_stress, pressure, temperature, dislocation_creep_parameters[i]);
                log_edot_and_deriv[i][1].first += log_volume_fraction;
              }
            if (use_peierls_creep)
              {
                log_edot_and_deriv[i][2] = peierls_creep->compute_approximate_log_strain_rate_and_derivative(log_viscoplastic_stress, pressure, temperature, peierls_creep_parameters[i]);
                log_edot_and_deriv[i][2].first += log_volume_fraction;
              }
            if (use_drucker_prager)
              {
                log_edot_and_deriv[i][3] = drucker_prager->compute_log_strain_rate_and_derivative(log_viscoplastic_stress, pressure, drucker_prager_parameters[i]);
                log_edot_and_deriv[i][3].first += log_volume_fraction;
              }
          }

        // 3) Calculate the total log strain rate from the first estimates
        // for the component log strain rates.
        // This will generally be more than the input total strain rate because
        // the creep stress in Step 1 was been calculated assuming that
        // only one mechanism was active, whereas the strain rate
        // calculated in Step 2 allowed all the mechanisms to
        // accommodate strain at that creep stress.
        std::pair<double, double> log_edot_ii_and_deriv_iterate = calculate_isostress_log_strain_rate_and_derivative(log_edot_and_deriv,
                                                                  viscoplastic_stress,
                                                                  inverse_kelvin_viscosity,
                                                                  partial_strain_rates);
        double log_strain_rate_residual = log_edot_ii_and_deriv_iterate.first - log_edot_ii;

        // 4) In this rheology model, the total strain rate is partitioned between
        // different flow components. We do not know how the strain is partitioned
        // between these components.

        // The following while loop contains a Newton iteration to obtain the
        // viscoplastic stress that is consistent with the total strain rate
        // and results in the correct partitioning of strain rate amongst
        // the different flow components.
        unsigned int stress_iteration = 0;
        while (std::abs(log_strain_rate_residual) > log_strain_rate_residual_threshold && stress_iteration < stress_max_iteration_number)
          {
            // Apply the Newton update for the log creep stress using the
            // strain-rate residual and strain-rate stress derivative
            double delta_log_viscoplastic_stress = log_strain_rate_residual / log_edot_ii_and_deriv_iterate.second;
            log_viscoplastic_stress += delta_log_viscoplastic_stress;
            viscoplastic_stress = std::exp(log_viscoplastic_stress);

            // Update the strain rates of all mechanisms with the new stress
            for (auto &log_edot_and_deriv_c : log_edot_and_deriv)
              {
                for (auto &i : active_flow_mechanisms)
                  log_edot_and_deriv_c[i].first += log_edot_and_deriv_c[i].second * delta_log_viscoplastic_stress;
              }

            // Compute the new log strain rate residual and log stress derivative
            log_edot_ii_and_deriv_iterate = calculate_isostress_log_strain_rate_and_derivative(log_edot_and_deriv,
                                            viscoplastic_stress,
                                            inverse_kelvin_viscosity,
                                            partial_strain_rates);
            log_strain_rate_residual = log_edot_ii - log_edot_ii_and_deriv_iterate.first;

            ++stress_iteration;

            // If anything that would be used in the next iteration is not finite, the
            // Newton iteration would trigger an exception and we want to abort the
            // iteration instead.
            // Currently, we still throw an exception, but if this exception is thrown,
            // another more robust iterative scheme should be implemented
            // (similar to that seen in the diffusion-dislocation material model).
            const bool abort_newton_iteration = !numbers::is_finite(log_viscoplastic_stress) || !numbers::is_finite(log_strain_rate_residual) || stress_iteration == stress_max_iteration_number;
            AssertThrow(!abort_newton_iteration,
                        ExcMessage("No convergence has been reached in the loop that determines "
                                   "the composite viscous creep stress. Aborting! "
                                   "Residual is " +
                                   Utilities::to_string(log_strain_rate_residual) +
                                   " after " + Utilities::to_string(stress_iteration) + " iterations. "
                                   "You can increase the number of iterations by adapting the "
                                   "parameter 'Maximum creep strain rate iterations'."));
          }

        // 5) We have now calculated the viscoplastic stress consistent with the total
        // strain rate, but the viscoplastic stress is only one component of the total stress,
        // because this material model also includes a viscosity damper
        // arranged in parallel with the viscoplastic elements.
        // The total stress is equal to the sum of the viscoplastic stress and
        // damper stress.
        double viscoplastic_strain_rate = 0.;
        for (auto &i : active_flow_mechanisms)
          viscoplastic_strain_rate += partial_strain_rates[i];

        const double damper_stress = 2. * damper_viscosity * viscoplastic_strain_rate;
        const double total_stress = viscoplastic_stress + damper_stress;

        // 6) Return the effective creep viscosity using the total stress
        return total_stress / (2. * edot_ii);
      }



      template <int dim>
      std::pair<double, double>
      CompositeViscoPlastic<dim>::calculate_isostress_log_strain_rate_and_derivative(const std::vector<std::array<std::pair<double, double>, 4>> &logarithmic_strain_rates_and_stress_derivatives,
                                                                                     const double viscoplastic_stress,
                                                                                     const double inverse_kelvin_viscosity,
                                                                                     std::vector<double> &partial_strain_rates) const
      {
        // The total strain rate
        double viscoplastic_strain_rate_sum = 0.0;

        // The sum of the stress derivatives multiplied by the mechanism strain rates
        double weighted_stress_derivative_sum = 0.0;

        // The first derivative of log(strain rate) with respect to log(stress)
        // is computed as sum_i(stress_exponent_i * edot_i) / sum_i(edot_i)
        // i.e., the stress exponents weighted by strain rate fraction
        // summed over the individual flow mechanisms (i).
        // Loop over active flow laws and add their contributions
        // to the strain rate and stress derivative

        // First, make sure that the active partial_strain_rates are reset to zero.
        for (auto &i : active_flow_mechanisms)
          partial_strain_rates[i] = 0;

        for (auto &logarithmic_strain_rates_and_stress_derivatives_c : logarithmic_strain_rates_and_stress_derivatives)
          {
            for (auto &i : active_flow_mechanisms)
              {
                double mechanism_log_strain_rate = logarithmic_strain_rates_and_stress_derivatives_c[i].first;

                // Check if the mechanism strain rate is within bounds to prevent underflow
                if (mechanism_log_strain_rate >= logmin)
                  {
                    const double mechanism_strain_rate = std::exp(mechanism_log_strain_rate);
                    partial_strain_rates[i] += mechanism_strain_rate;
                    const double log_stress_derivative = logarithmic_strain_rates_and_stress_derivatives_c[i].second;
                    viscoplastic_strain_rate_sum += mechanism_strain_rate;
                    weighted_stress_derivative_sum += log_stress_derivative * mechanism_strain_rate;
                  }
              }
          }

        const double log_viscoplastic_strain_rate_derivative = weighted_stress_derivative_sum / viscoplastic_strain_rate_sum;

        // Some opaque mathematics converts the viscoplastic strain rate to the total strain rate.
        const double inverse_hard_and_elastic_viscosity = (inverse_maximum_viscosity + inverse_kelvin_viscosity);
        const double f = 0.5 * inverse_hard_and_elastic_viscosity * viscoplastic_stress;
        const double g = 1. + minimum_viscosity * inverse_kelvin_viscosity;
        const double strain_rate = (strain_rate_scaling_factor * g * viscoplastic_strain_rate_sum) + f;
        const double strain_rate_hard_and_elastic = strain_rate - viscoplastic_strain_rate_sum;

        // Elastic strain rate
        partial_strain_rates[elastic_strain_rate_index] = strain_rate_hard_and_elastic * (inverse_kelvin_viscosity / inverse_hard_and_elastic_viscosity);

        // Hard damper strain rate
        partial_strain_rates[damper_strain_rate_index] = strain_rate_hard_and_elastic - partial_strain_rates[elastic_strain_rate_index];

        // And the partial derivative of the log *total* strain rate
        // with respect to log *viscoplastic* stress follows as
        const double log_strain_rate_derivative = (strain_rate_scaling_factor * g * viscoplastic_strain_rate_sum * log_viscoplastic_strain_rate_derivative + f) / strain_rate;

        return std::make_pair(std::log(strain_rate), log_strain_rate_derivative);
      }



      template <int dim>
      double
      CompositeViscoPlastic<dim>::compute_isostrain_viscosity (const double pressure,
                                                               const double temperature,
                                                               const double grain_size,
                                                               const std::vector<double> &volume_fractions,
                                                               const SymmetricTensor<2,dim> &strain_rate,
                                                               std::vector<double> &partial_strain_rates,
                                                               const std::vector<double> &phase_function_values,
                                                               const std::vector<unsigned int> &n_phase_transitions_per_composition) const
      {
        // The isostrain viscosity is calculated by summing the viscosities
        //of each composition weighted by their volume fractions.
        double viscosity = 0.;
        double total_volume_fraction = 1.;
        for (unsigned int composition=0; composition < number_of_chemical_compositions; ++composition)
          {
            // Only include the contribution to the viscosity
            // from a given composition if the volume fraction exceeds
            // a certain (small) fraction.
            if (volume_fractions[composition] > 2. * std::numeric_limits<double>::epsilon())
              {
                // There is no elastic component allowed in isostrain models,
                // so the size of vector partial_strain_rates_composition is
                // one smaller than partial_strain_rates.
                std::vector<double> partial_strain_rates_composition(n_decomposed_strain_rates-1, 0.);
                viscosity += (volume_fractions[composition]
                              * compute_composition_viscosity (pressure,
                                                               temperature,
                                                               grain_size,
                                                               composition,
                                                               strain_rate,
                                                               partial_strain_rates_composition,
                                                               phase_function_values,
                                                               n_phase_transitions_per_composition));
                for (auto &j : active_flow_mechanisms)
                  partial_strain_rates[j] += volume_fractions[composition] * partial_strain_rates_composition[j];

                // Shift the strain rate for the hard viscosity damper into the last position in partial_strain_rates.
                partial_strain_rates[damper_strain_rate_index] += volume_fractions[composition] * partial_strain_rates_composition[isostrain_damper_strain_rate_index];
              }
            else
              {
                total_volume_fraction -= volume_fractions[composition];
              }
          }

        viscosity /= total_volume_fraction;

        for (auto &j : active_flow_mechanisms)
          partial_strain_rates[j] /= total_volume_fraction;

        return viscosity;
      }



      template <int dim>
      double
      CompositeViscoPlastic<dim>::compute_composition_viscosity (const double pressure,
                                                                 const double temperature,
                                                                 const double grain_size,
                                                                 const unsigned int composition,
                                                                 const SymmetricTensor<2,dim> &strain_rate,
                                                                 std::vector<double> &partial_strain_rates,
                                                                 const std::vector<double> &phase_function_values,
                                                                 const std::vector<unsigned int> &n_phase_transitions_per_composition) const
      {
        // Calculate the viscosity for a single compositional field based on the local state.
        // If strain rate is zero (like during the first time step) set it to some very small number
        // to prevent a division-by-zero, and a floating point exception.
        // Otherwise, calculate the square-root of the norm of the second invariant of the deviatoric-
        // strain rate (often simplified as epsilondot_ii)
        const double edot_ii = std::max(std::sqrt(std::max(-second_invariant(deviator(strain_rate)), 0.)),
                                        minimum_strain_rate);
        const double log_edot_ii = std::log(edot_ii);

        Rheology::DiffusionCreepParameters diffusion_creep_parameters;
        Rheology::DislocationCreepParameters dislocation_creep_parameters;
        Rheology::PeierlsCreepParameters peierls_creep_parameters;
        Rheology::DruckerPragerParameters drucker_prager_parameters;

        // 1) Estimate the stress running through the creep elements and the
        // maximum viscosity element, whose viscosity is defined in parse_parameters.
        // These are all arranged in series. Taking the minimum viscosity from
        // all of these elements provides an excellent first approximation
        // to the true viscosity.
        // The stress can then be calculated as 2 * viscoplastic_viscosity_guess * edot_ii
        double viscoplastic_viscosity_guess = maximum_viscosity;

        if (use_diffusion_creep)
          {
            diffusion_creep_parameters = diffusion_creep->compute_creep_parameters(composition, phase_function_values, n_phase_transitions_per_composition);
            viscoplastic_viscosity_guess = std::min(diffusion_creep->compute_viscosity(pressure, temperature, grain_size, composition, phase_function_values, n_phase_transitions_per_composition), viscoplastic_viscosity_guess);
          }

        if (use_dislocation_creep)
          {
            dislocation_creep_parameters = dislocation_creep->compute_creep_parameters(composition, phase_function_values, n_phase_transitions_per_composition);
            viscoplastic_viscosity_guess = std::min(dislocation_creep->compute_viscosity(edot_ii, pressure, temperature, composition, phase_function_values, n_phase_transitions_per_composition), viscoplastic_viscosity_guess);
          }

        if (use_peierls_creep)
          {
            peierls_creep_parameters = peierls_creep->compute_creep_parameters(composition);
            viscoplastic_viscosity_guess = std::min(peierls_creep->compute_approximate_viscosity(edot_ii, pressure, temperature, composition), viscoplastic_viscosity_guess);
          }

        if (use_drucker_prager)
          {
            drucker_prager_parameters = drucker_prager->compute_drucker_prager_parameters(composition, phase_function_values, n_phase_transitions_per_composition);
            viscoplastic_viscosity_guess = std::min(drucker_prager->compute_viscosity(drucker_prager_parameters.cohesion,
                                                                                      drucker_prager_parameters.angle_internal_friction, pressure, edot_ii, drucker_prager_parameters.max_yield_stress), viscoplastic_viscosity_guess);
          }

        double viscoplastic_stress = 2. * viscoplastic_viscosity_guess * edot_ii;
        double log_viscoplastic_stress = std::log(viscoplastic_stress);

        // 2) Calculate the log strain rate and first derivative with respect to
        // log stress for each element using the guessed creep stress.
        std::array<std::pair<double, double>, 4> log_edot_and_deriv;

        if (use_diffusion_creep)
          log_edot_and_deriv[0] = diffusion_creep->compute_log_strain_rate_and_derivative(log_viscoplastic_stress, pressure, temperature, grain_size, diffusion_creep_parameters);

        if (use_dislocation_creep)
          log_edot_and_deriv[1] = dislocation_creep->compute_log_strain_rate_and_derivative (log_viscoplastic_stress, pressure, temperature, dislocation_creep_parameters);

        if (use_peierls_creep)
          log_edot_and_deriv[2] = peierls_creep->compute_approximate_log_strain_rate_and_derivative(log_viscoplastic_stress, pressure, temperature, peierls_creep_parameters);

        if (use_drucker_prager)
          log_edot_and_deriv[3] = drucker_prager->compute_log_strain_rate_and_derivative (log_viscoplastic_stress, pressure, drucker_prager_parameters);

        // 3) Calculate the total log strain rate from the first estimates
        // for the component log strain rates.
        // This will generally be more than the input total strain rate because
        // the creep stress in Step 1 was been calculated assuming that
        // only one mechanism was active, whereas the strain rate
        // calculated in Step 2 allowed all the mechanisms to
        // accommodate strain at that creep stress.
        std::pair<double, double> log_edot_ii_and_deriv_iterate = calculate_composition_log_strain_rate_and_derivative(log_edot_and_deriv,
                                                                  viscoplastic_stress,
                                                                  partial_strain_rates);
        double log_strain_rate_residual = log_edot_ii_and_deriv_iterate.first - log_edot_ii;

        // 4) In this rheology model, the total strain rate is partitioned between
        // different flow components. We do not know how the strain is partitioned
        // between these components.

        // The following while loop contains a Newton iteration to obtain the
        // viscoplastic stress that is consistent with the total strain rate.
        unsigned int stress_iteration = 0;
        while (std::abs(log_strain_rate_residual) > log_strain_rate_residual_threshold
               && stress_iteration < stress_max_iteration_number)
          {
            // Apply the Newton update for the log creep stress using the
            // strain-rate residual and strain-rate stress derivative
            double delta_log_viscoplastic_stress = log_strain_rate_residual/log_edot_ii_and_deriv_iterate.second;
            log_viscoplastic_stress += delta_log_viscoplastic_stress;
            viscoplastic_stress = std::exp(log_viscoplastic_stress);

            // Update the strain rates of all mechanisms with the new stress
            for (auto &i : active_flow_mechanisms)
              log_edot_and_deriv[i].first += log_edot_and_deriv[i].second * delta_log_viscoplastic_stress;

            // Compute the new log strain rate residual and log stress derivative
            log_edot_ii_and_deriv_iterate = calculate_composition_log_strain_rate_and_derivative(log_edot_and_deriv,
                                            viscoplastic_stress,
                                            partial_strain_rates);
            log_strain_rate_residual = log_edot_ii - log_edot_ii_and_deriv_iterate.first;

            ++stress_iteration;

            // If anything that would be used in the next iteration is not finite, the
            // Newton iteration would trigger an exception and we want to abort the
            // iteration instead.
            // Currently, we still throw an exception, but if this exception is thrown,
            // another more robust iterative scheme should be implemented
            // (similar to that seen in the diffusion-dislocation material model).
            const bool abort_newton_iteration = !numbers::is_finite(log_viscoplastic_stress)
                                                || !numbers::is_finite(log_strain_rate_residual)
                                                || stress_iteration == stress_max_iteration_number;
            AssertThrow(!abort_newton_iteration,
                        ExcMessage("No convergence has been reached in the loop that determines "
                                   "the composite viscous creep stress. Aborting! "
                                   "Residual is " + Utilities::to_string(log_strain_rate_residual) +
                                   " after " + Utilities::to_string(stress_iteration) + " iterations. "
                                   "You can increase the number of iterations by adapting the "
                                   "parameter 'Maximum creep strain rate iterations'."));
          }

        // 5) We have now calculated the viscoplastic stress consistent with the total
        // strain rate, but the viscoplastic stress is only one component of the total stress,
        // because this material model also includes a viscosity damper
        // arranged in parallel with the viscoplastic elements.
        // The total stress is equal to the sum of the viscoplastic stress and
        // damper stress.
        const double damper_stress = 2. * damper_viscosity * (edot_ii - partial_strain_rates[isostrain_damper_strain_rate_index]);
        const double total_stress = viscoplastic_stress + damper_stress;

        // 6) Return the effective creep viscosity using the total stress
        return total_stress/(2.*edot_ii);
      }



      template <int dim>
      std::pair<double, double>
      CompositeViscoPlastic<dim>::calculate_composition_log_strain_rate_and_derivative(const std::array<std::pair<double, double>, 4> &logarithmic_strain_rates_and_stress_derivatives,
          const double viscoplastic_stress,
          std::vector<double> &partial_strain_rates) const
      {
        // Initialize the double containing the total strain rate
        // for the single compositional field of interest.
        double viscoplastic_strain_rate_sum = 0.0;

        // The sum of the stress derivatives multiplied by the mechanism strain rates
        double weighted_stress_derivative_sum = 0.0;

        // The first derivative of log(strain rate) with respect to log(stress)
        // is computed as sum_i(stress_exponent_i * edot_i) / sum_i(edot_i)
        // i.e., the stress exponents weighted by strain rate fraction
        // summed over the individual flow mechanisms (i).

        // Loop over active flow laws and add their contributions
        // to the strain rate and stress derivative
        for (auto &i : active_flow_mechanisms)
          {
            double mechanism_log_strain_rate = logarithmic_strain_rates_and_stress_derivatives[i].first;

            // Check if the mechanism strain rate is within bounds to prevent underflow
            if (mechanism_log_strain_rate >= logmin)
              {
                const double mechanism_strain_rate = std::exp(mechanism_log_strain_rate);
                partial_strain_rates[i] = mechanism_strain_rate;
                const double log_stress_derivative = logarithmic_strain_rates_and_stress_derivatives[i].second;
                viscoplastic_strain_rate_sum += mechanism_strain_rate;
                weighted_stress_derivative_sum += log_stress_derivative * mechanism_strain_rate;
              }
          }

        const double log_viscoplastic_strain_rate_derivative = weighted_stress_derivative_sum / viscoplastic_strain_rate_sum;

        // Some opaque mathematics converts the viscoplastic strain rate to the total strain rate.
        const double f = viscoplastic_stress / (2. * maximum_viscosity);
        const double strain_rate = (strain_rate_scaling_factor * viscoplastic_strain_rate_sum) + f;
        partial_strain_rates[isostrain_damper_strain_rate_index] = strain_rate - viscoplastic_strain_rate_sum;
        // And the partial derivative of the log *total* strain rate
        // with respect to log *viscoplastic* stress follows as
        const double log_strain_rate_derivative = (strain_rate_scaling_factor * viscoplastic_strain_rate_sum * log_viscoplastic_strain_rate_derivative + f) / strain_rate;

        return std::make_pair(std::log(strain_rate), log_strain_rate_derivative);
      }



      // TODO: UNCOMMENT
      /*
      template <int dim>
      void
      CompositeViscoPlastic<dim>::create_elastic_additional_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        // Create the ElasticAdditionalOutputs that include the average shear modulus, elastic
        // viscosity, timestep ratio and total deviatoric stress of the current timestep.
        if (out.template get_additional_output<ElasticAdditionalOutputs<dim>>() == nullptr)
          {
            const unsigned int n_points = out.n_evaluation_points();
            out.additional_outputs.push_back(
              std::make_unique<ElasticAdditionalOutputs<dim>> (n_points));
          }

        // We need to modify the shear heating outputs to correctly account for elastic stresses.
        if (out.template get_additional_output<HeatingModel::PrescribedShearHeatingOutputs<dim>>() == nullptr)
          {
            const unsigned int n_points = out.n_evaluation_points();
            out.additional_outputs.push_back(
              std::make_unique<HeatingModel::PrescribedShearHeatingOutputs<dim>> (n_points));
          }

        // Create the ReactionRateOutputs that are necessary for the operator splitting
        // step (either on the fields or directly on the particles)
        // that sets both sets of stresses to the total stress of the
        // previous timestep.
        if (out.template get_additional_output<ReactionRateOutputs<dim>>() == nullptr &&
            (this->get_parameters().use_operator_splitting || (this->get_parameters().mapped_particle_properties).count(this->introspection().compositional_index_for_name("ve_stress_xx"))))
          {
            const unsigned int n_points = out.n_evaluation_points();
            out.additional_outputs.push_back(
              std::make_unique<MaterialModel::ReactionRateOutputs<dim>>(n_points, this->n_compositional_fields()));
          }
      }



      template <int dim>
      void
      CompositeViscoPlastic<dim>::fill_elastic_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
                                                        const std::vector<double> &inverse_kelvin_viscosities,
                                                        MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        // Create a reference to the structure for the elastic outputs.
        // The structure is created during the Stokes assembly.
        MaterialModel::ElasticOutputs<dim>
        *elastic_out = out.template get_additional_output<MaterialModel::ElasticOutputs<dim>>();

        // Create a reference to the structure for the prescribed shear heating outputs.
        // The structure is created during the advection assembly.
        HeatingModel::PrescribedShearHeatingOutputs<dim>
        *heating_out = out.template get_additional_output<HeatingModel::PrescribedShearHeatingOutputs<dim>>();

        if (elastic_out == nullptr && heating_out == nullptr)
          return;

        // TODO should a RHS term be a separate MaterialProperties?
        if (in.requests_property(MaterialProperties::additional_outputs))
          {
            // The viscosity should be averaged if material averaging is applied.
            std::vector<double> effective_creep_viscosities;
            if (this->get_parameters().material_averaging != MaterialAveraging::none)
              {
                MaterialModelOutputs<dim> out_copy(out.n_evaluation_points(),
                                                   this->introspection().n_compositional_fields);
                out_copy.viscosities = out.viscosities;

                const MaterialAveraging::AveragingOperation averaging_operation_for_viscosity =
                  get_averaging_operation_for_viscosity(this->get_parameters().material_averaging);
                MaterialAveraging::average(averaging_operation_for_viscosity,
                                           in.current_cell,
                                           this->introspection().quadratures.velocities,
                                           this->get_mapping(),
                                           in.requested_properties,
                                           out_copy);

                effective_creep_viscosities = out_copy.viscosities;
              }
            else
              effective_creep_viscosities = out.viscosities;

            const unsigned int stress_start_index = this->introspection().compositional_index_for_name("ve_stress_xx");

            for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
              {
                const SymmetricTensor<2, dim> deviatoric_strain_rate = deviator(in.strain_rate[i]);

                // Get stress from timestep $t$ rotated and advected into the current
                // timestep $t+\Delta t_c$ from the compositional fields.
                // This function is only evaluated during the assembly of the Stokes equations
                // (the force term goes into the rhs of the momentum equation).
                // This happens after the advection equations have been solved, and hence in.composition
                // contains the rotated and advected stresses $tau^{0adv}$.
                // Only at the beginning of the next timestep do we add the stress update of the
                // current timestep to the stress stored in the compositional fields, giving
                // $\tau{t+\Delta t_c}$ with $t+\Delta t_c$ being the current timestep.
                const SymmetricTensor<2,dim> stress_0_advected (Utilities::Tensors::to_symmetric_tensor<dim>(&in.composition[i][stress_start_index],
                                                                &in.composition[i][stress_start_index]+n_independent_components));

                // Average effective creep viscosity
                // Use the viscosity corresponding to the stresses selected above.
                // out.viscosities is computed during the assembly of the Stokes equations
                // based on the current_linearization_point. This means that it will be updated after every
                // nonlinear Stokes iteration.
                // The effective creep viscosity has already been scaled with the timestep ratio dtc/dte.
                const double effective_creep_viscosity = effective_creep_viscosities[i];

                // The force term is computed as:
                // $\frac{-\eta_{effcreep} \tau_{0adv}}{\eta_{e}}$, where $\eta_{effcreep}$ is the
                // current harmonic average of the viscous and elastic viscosity, or the yield stress
                // divided by two times the second invariant of the deviatoric strain rate.
                // In case the computational timestep differs from the elastic timestep,
                // linearly interpolate between the two.
                // The elastic viscosity has also already been scaled with the timestep ratio.
                const double viscosity_ratio = effective_creep_viscosity * inverse_kelvin_viscosities[i];

                if (elastic_out != nullptr)
                  {
                    elastic_out->elastic_force[i] = -1. * viscosity_ratio * stress_0_advected;

                    // The viscoelastic strain rate is needed only when the Newton method is selected.
                    const typename Parameters<dim>::NonlinearSolver::Kind nonlinear_solver = this->get_parameters().nonlinear_solver;
                    if ((nonlinear_solver == Parameters<dim>::NonlinearSolver::iterated_Advection_and_Newton_Stokes) ||
                        (nonlinear_solver == Parameters<dim>::NonlinearSolver::single_Advection_iterated_Newton_Stokes))
                      elastic_out->viscoelastic_strain_rate[i] = compute_effective_strain_rate(in.strain_rate[i], stress_0_advected, inverse_kelvin_viscosities[i]);
                  }

                // The shear heating rate (used by the heating model) depends not only
                // on the stress and strain rates from the viscous flow laws, but also from
                // the three viscous dampers. Usually, the heating from the dampers will be
                // minor, but for consistency they should be included.

                // 1) The total stress
                const SymmetricTensor<2, dim> stress = 2. * effective_creep_viscosity * deviatoric_strain_rate + viscosity_ratio * stress_0_advected;

                // 2) The elastic damper
                const SymmetricTensor<2, dim> kelvin_strain_rate = 0.5 * (stress - stress_0_advected) * inverse_kelvin_viscosities[i];
                const double damper_viscosity = elasticity->get_damper_viscosity();
                const double damper_power_density = 2. * damper_viscosity * kelvin_strain_rate * kelvin_strain_rate;

                // 3) Total other work (assuming incompressibility)
                // This includes the power density from the other two dampers
                const SymmetricTensor<2, dim> visco_plastic_strain_rate = deviatoric_strain_rate - kelvin_strain_rate;
                const double viscoplastic_power_density = stress * visco_plastic_strain_rate;
                // If compressible,
                // visco_plastic_strain_rate = visco_plastic_strain_rate -
                //                             1. / 3. * trace(visco_plastic_strain_rate) * unit_symmetric_tensor<dim>();

                // The shear heating term needs to account for the elastic stress, but only the visco_plastic strain rate.
                // This is best computed here, and stored for later use by the heating model.
                if (heating_out != nullptr)
                  heating_out->prescribed_shear_heating_rates[i] = viscoplastic_power_density + damper_power_density;
              }

          }
      }



      template <int dim>
      void
      CompositeViscoPlastic<dim>::fill_elastic_additional_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
                                                                   const std::vector<double> &inverse_kelvin_viscosities,
                                                                   MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        // Create a reference to the structure for the elastic additional outputs
        MaterialModel::ElasticAdditionalOutputs<dim>
        *elastic_additional_out = out.template get_additional_output<MaterialModel::ElasticAdditionalOutputs<dim>>();

        if (elastic_additional_out == nullptr || !in.requests_property(MaterialProperties::additional_outputs))
          return;

        const unsigned int stress_start_index = this->introspection().compositional_index_for_name("ve_stress_xx");
        const double dtc = this->get_timestep();
        const double elastic_damper_viscosity = elasticity->get_damper_viscosity();

        for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
          {
            const double effective_viscosity = out.viscosities[i];
            const SymmetricTensor<2, dim> deviatoric_strain_rate = deviator(in.strain_rate[i]);
            const SymmetricTensor<2,dim> stress_0_advected (Utilities::Tensors::to_symmetric_tensor<dim>(&in.composition[i][stress_start_index],
                                                            &in.composition[i][stress_start_index]+n_independent_components));

            // Apply the stress update to get the total deviatoric stress of timestep t.
            // This is not the stress passing through the elastic element, because of the elastic damper.
            elastic_additional_out->deviatoric_stress[i] = 2. * effective_viscosity * deviatoric_strain_rate + effective_viscosity * inverse_kelvin_viscosities[i] * stress_0_advected;
            elastic_additional_out->elastic_viscosity[i] = 1. / inverse_kelvin_viscosities[i];
            elastic_additional_out->elastic_shear_moduli[i] = (elastic_additional_out->elastic_viscosity[i] - elastic_damper_viscosity)/dtc;
          }
      }
      */


      // Rotate the elastic stresses of the previous timestep $t$ into the current timestep $t+dtc$.
      template <int dim>
      void
      CompositeViscoPlastic<dim>::fill_reaction_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
                                                         const std::vector<double> &,
                                                         MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        elasticity->fill_reaction_outputs(in, std::vector<double>(), out);
      }



      // The following function computes the reaction rates for the operator
      // splitting step that at the beginning of the new timestep $t+dtc$ updates the
      // stored compositions $tau^{0\mathrm{adv}}$ at time $t$ to $tau^{t}$.
      // This update consists of the stress change resulting from system evolution,
      // but does not advect or rotate the elastic stress tensor. Advection is done by
      // solving the advection equation and the elastic stress tensor is rotated through
      // the source term (reaction_terms) of that same equation.
      template <int dim>
      void
      CompositeViscoPlastic<dim>::fill_reaction_rates (const MaterialModel::MaterialModelInputs<dim> &in,
                                                       const std::vector<double> &inverse_kelvin_viscosities,
                                                       MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        ReactionRateOutputs<dim> *reaction_rate_out = out.template get_additional_output<ReactionRateOutputs<dim>>();

        if (reaction_rate_out == nullptr)
          return;

        // At the moment when the reaction rates are required (at the beginning of the timestep),
        // the solution vector 'solution' holds the stress from the previous timestep,
        // advected into the new position of the previous timestep, so $\tau^{t}_{0adv}$.
        // This is the same as the vector 'old_solution' holds. At later moments during the current timestep,
        // 'solution' will hold the current_linearization_point instead of the solution of the previous timestep.
        //
        // In case fields are used to track the stresses, MaterialModelInputs are based on 'solution'
        // when calling the MaterialModel for the reaction rates. When particles are used, MaterialModelInputs
        // for this function are filled with the old_solution (including for the strain rate), except for the
        // compositions that represent the stress tensor components, these are taken directly from the
        // particles. As the particles are restored to their pre-advection location at the beginning of
        // each nonlinear iteration, their values and positions correspond to the old solution.
        // This means that in both cases we can use 'in' to get to the $\tau^{t}_{0adv}$ and velocity/strain rate of the
        // previous timestep.
        if (in.current_cell.state() == IteratorState::valid && this->get_timestep_number() > 0 && in.requests_property(MaterialProperties::reaction_rates))
          {
            const unsigned int stress_start_index = this->introspection().compositional_index_for_name("ve_stress_xx");
            const double elastic_damper_viscosity = elasticity->get_damper_viscosity();
            for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
              {
                // Set all reaction rates to zero
                for (unsigned int c = 0; c < in.composition[i].size(); ++c)
                  reaction_rate_out->reaction_rates[i][c] = 0.0;

                // Get $\tau^{0adv}$ of the previous timestep t from the compositional fields.
                // This stress includes the rotation and advection of the previous timestep,
                // i.e., the reaction term (which prescribes the change in stress due to rotation
                // over the previous timestep) has already been applied during the previous timestep.
                const SymmetricTensor<2, dim> stress_0_t (Utilities::Tensors::to_symmetric_tensor<dim>(&in.composition[i][stress_start_index],
                                                          &in.composition[i][stress_start_index]+n_independent_components));

                // $\eta^{t}_{effcreep}$. This viscosity has been calculated with the timestep_ratio dtc/dte.
                const double effective_creep_viscosity = out.viscosities[i];

                // Compute the total stress at time t.
                const SymmetricTensor<2, dim>
                stress_t = 2. * effective_creep_viscosity * deviator(in.strain_rate[i])
                           + effective_creep_viscosity * inverse_kelvin_viscosities[i] * stress_0_t;

                // Fill reaction rates.
                // During this timestep, the reaction rates will be multiplied
                // with the current timestep size to turn the rate of change into a change.
                // However, this update belongs
                // to the previous timestep. Therefore we divide by the
                // current timestep and multiply with the previous one.
                // When multiplied with the current timestep, this will give
                // (rate * previous_dt / current_dt) * current_dt = rate * previous_dt = previous_change.
                // previous_change = (1 - eta_d/eta_kel)*(stress_t - stress_0_t).
                // To compute the rate we should return to the operator splitting scheme,
                // we therefore divide the change in stress by the current timestep current_dt (= dtc).

                const SymmetricTensor<2, dim> stress_update = (1. - (elastic_damper_viscosity*inverse_kelvin_viscosities[i])) * (stress_t - stress_0_t) / this->get_timestep();

                Utilities::Tensors::unroll_symmetric_tensor_into_array(stress_update,
                                                                       &reaction_rate_out->reaction_rates[i][stress_start_index],
                                                                       &reaction_rate_out->reaction_rates[i][stress_start_index]+n_independent_components);
              }
          }
      }



      // Overload the + operator to act on two pairs of doubles.
      std::pair<double,double> operator+(const std::pair<double,double> &x, const std::pair<double,double> &y)
      {
        return std::make_pair(x.first+y.first, x.second+y.second);
      }



      template <int dim>
      void
      CompositeViscoPlastic<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Viscosity averaging scheme", "isostress",
                           Patterns::Selection("isostress|isostrain"),
                           "Determines the relationship between the conditions experienced by "
                           "each chemical compositional field in a composite material. "
                           "Select either isostress (the default option, Reuss averaging) "
                           "or isostrain (Voigt averaging).");

        prm.declare_entry ("Include diffusion creep in composite rheology", "true",
                           Patterns::Bool (),
                           "Whether to include diffusion creep in the composite rheology formulation.");

        prm.declare_entry ("Include dislocation creep in composite rheology", "true",
                           Patterns::Bool (),
                           "Whether to include dislocation creep in the composite rheology formulation.");

        prm.declare_entry ("Include Peierls creep in composite rheology", "true",
                           Patterns::Bool (),
                           "Whether to include Peierls creep in the composite rheology formulation.");

        prm.declare_entry ("Include Drucker Prager plasticity in composite rheology", "true",
                           Patterns::Bool (),
                           "Whether to include Drucker-Prager plasticity in the composite rheology formulation.");

        prm.declare_entry ("Include elasticity in composite rheology", "false",
                           Patterns::Bool (),
                           "Whether to include elasticity in the composite rheology formulation.");

        // Diffusion creep parameters
        Rheology::DiffusionCreep<dim>::declare_parameters(prm);

        // Dislocation creep parameters
        Rheology::DislocationCreep<dim>::declare_parameters(prm);

        // Dislocation creep parameters
        Rheology::PeierlsCreep<dim>::declare_parameters(prm);

        // Drucker Prager parameters
        Rheology::DruckerPragerPower<dim>::declare_parameters(prm);

        // Elastic parameters
        Rheology::Elasticity<dim>::declare_parameters(prm);

        // Some of the parameters below are shared with the subordinate
        // rheology models (diffusion, dislocation, ...),
        // and will already have been declared. This is fine, the deal.II
        // parameter handler allows repeated declarations. The default values
        // below will override any defaults given in the subordinate
        // rheology models. The new defaults will apply both to
        // this rheology model and to the subordinate rheology modules.
        prm.declare_entry ("Minimum strain rate", "1.4e-20", Patterns::Double(0.),
                           "Stabilizes strain dependent viscosity. Units: \\si{\\per\\second}.");

        // Viscosity iteration parameters
        prm.declare_entry ("Strain rate residual tolerance", "1e-10", Patterns::Double(0.),
                           "Tolerance for correct log strain rate residual.");
        prm.declare_entry ("Maximum creep strain rate iterations", "40", Patterns::Integer(0),
                           "Maximum number of iterations to find the correct "
                           "viscoplastic strain rate.");

        // Strain rate and stress limiting parameters
        prm.declare_entry ("Minimum viscosity", "1.e17",
                           Patterns::Double(0.),
                           "Minimum effective viscosity. Units: \\si{\\pascal\\second}.");

        prm.declare_entry ("Maximum viscosity", "1.e28",
                           Patterns::Double(0.),
                           "Maximum effective viscosity. Units: \\si{\\pascal\\second}.");
      }



      template <int dim>
      void
      CompositeViscoPlastic<dim>::parse_parameters (ParameterHandler &prm,
                                                    const std::unique_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition)
      {
        viscosity_averaging_scheme = ViscosityAveraging::parse("Viscosity averaging scheme",
                                                               prm);

        // A background field is required by the subordinate material models
        number_of_chemical_compositions = this->introspection().n_chemical_composition_fields() + 1;

        minimum_strain_rate = prm.get_double("Minimum strain rate");

        // Iteration parameters
        log_strain_rate_residual_threshold = prm.get_double ("Strain rate residual tolerance");
        stress_max_iteration_number = prm.get_integer ("Maximum creep strain rate iterations");

        // Read maximum viscosity parameter
        maximum_viscosity = prm.get_double ("Maximum viscosity");
        inverse_maximum_viscosity = 1. / maximum_viscosity;

        // Process minimum viscosity parameter
        // In this rheology model, there are two viscous dampers designed
        // to stabilise the solution, a "limiting damper" arranged in parallel
        // with the flow law components which stops the viscosity going to zero,
        // and a "maximum viscosity damper", which is placed in series with the
        // combined flow law components and the limiting damper.
        // - If the real creep viscosity is infinite, the total strain rate will
        // run through only the "maximum viscosity" damper.
        // - If the real creep viscosity is equal to zero, the total strain rate
        // will run through the "limiting damper" and the "maximum viscosity" damper,
        // with a total viscosity equal to the harmonic sum of these dampers.
        // Therefore, the "limiting damper" has a viscosity equal to
        // eta_max * eta_min / (eta_max - eta_min).
        // When scaling the viscoplastic strain up to the total strain,
        // eta_max / (eta_max - eta_min) becomes a useful value,
        // which we here call the "strain_rate_scaling_factor".
        minimum_viscosity = prm.get_double("Minimum viscosity");

        AssertThrow(minimum_viscosity > 0,
                    ExcMessage("Minimum viscosity needs to be larger than zero."));

        AssertThrow(maximum_viscosity > 1.1 * minimum_viscosity,
                    ExcMessage("The maximum viscosity needs to be at least ten percent larger than the minimum viscosity. "
                               "If you require an isoviscous model consider a different rheology, or set the "
                               "parameters of the active flow laws to be independent of temperature, pressure, grain size, and stress."));

        strain_rate_scaling_factor = maximum_viscosity / (maximum_viscosity - minimum_viscosity);
        damper_viscosity = maximum_viscosity * minimum_viscosity / (maximum_viscosity - minimum_viscosity);

        // Rheological parameters
        // Diffusion creep parameters
        use_diffusion_creep = prm.get_bool("Include diffusion creep in composite rheology");
        if (use_diffusion_creep)
          {
            diffusion_creep = std::make_unique<Rheology::DiffusionCreep<dim>>();
            diffusion_creep->initialize_simulator (this->get_simulator());
            diffusion_creep->parse_parameters(prm, expected_n_phases_per_composition);

            active_flow_mechanisms.push_back(0);
          }

        // Dislocation creep parameters
        use_dislocation_creep = prm.get_bool ("Include dislocation creep in composite rheology");
        if (use_dislocation_creep)
          {
            dislocation_creep = std::make_unique<Rheology::DislocationCreep<dim>>();
            dislocation_creep->initialize_simulator (this->get_simulator());
            dislocation_creep->parse_parameters(prm, expected_n_phases_per_composition);

            active_flow_mechanisms.push_back(1);
          }

        // Peierls creep parameters
        use_peierls_creep = prm.get_bool ("Include Peierls creep in composite rheology");
        if (use_peierls_creep)
          {
            peierls_creep = std::make_unique<Rheology::PeierlsCreep<dim>>();
            peierls_creep->initialize_simulator (this->get_simulator());
            peierls_creep->parse_parameters(prm, expected_n_phases_per_composition);
            active_flow_mechanisms.push_back(2);

            AssertThrow((prm.get ("Peierls creep flow law") == "viscosity approximation"),
                        ExcMessage("The Peierls creep flow law parameter needs to be set to viscosity approximation."));

          }

        // Drucker Prager parameters
        use_drucker_prager = prm.get_bool ("Include Drucker Prager plasticity in composite rheology");
        if (use_drucker_prager)
          {
            drucker_prager = std::make_unique<Rheology::DruckerPragerPower<dim>>();
            drucker_prager->initialize_simulator (this->get_simulator());
            drucker_prager->parse_parameters(prm, expected_n_phases_per_composition);
            active_flow_mechanisms.push_back(3);
          }

        AssertThrow(active_flow_mechanisms.size() > 0,
                    ExcMessage("You need to include at least one non-elastic deformation mechanism."));

        // Elastic parameters
        // Elasticity is not treated as a flow mechanism,
        // so we do not push back the active_flow_mechanisms vector.
        use_elasticity = prm.get_bool ("Include elasticity in composite rheology");
        if (use_elasticity)
          {
            elasticity = std::make_unique<Rheology::Elasticity<dim>>();
            elasticity->initialize_simulator (this->get_simulator());
            elasticity->parse_parameters(prm);
            AssertThrow(viscosity_averaging_scheme == ViscosityAveraging::isostress,
                        ExcMessage("Elasticity in the CompositeViscoPlastic rheology "
                                   "requires that the 'Viscosity averaging scheme' be "
                                   "set to isostress."));
            AssertThrow(prm.get ("Use fixed elastic time step") == "false",
                        ExcMessage("Elasticity in the CompositeViscoPlastic rheology "
                                   "requires that 'Use fixed elastic time step' be "
                                   "set to false."));
            AssertThrow(std::abs(prm.get_double ("Stabilization time scale factor") - 1.) < std::numeric_limits<double>::min(),
                        ExcMessage("Elasticity in the CompositeViscoPlastic rheology "
                                   "requires that the 'Stabilization time scale factor' be "
                                   "set to 1."));
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
  namespace Rheology \
  { \
    template class CompositeViscoPlastic<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
