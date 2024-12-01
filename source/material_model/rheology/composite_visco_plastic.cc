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
      // Drucker-Prager plasticity and a constant (high) viscosity limiter.
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
      CompositeViscoPlastic<dim>::compute_viscosity (const double pressure,
                                                     const double temperature,
                                                     const double grain_size,
                                                     const std::vector<double> &volume_fractions,
                                                     const SymmetricTensor<2,dim> &strain_rate,
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
                                                       strain_rate,
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
                                                       strain_rate,
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
                                                               const SymmetricTensor<2,dim> &strain_rate,
                                                               std::vector<double> &partial_strain_rates,
                                                               const std::vector<double> &phase_function_values,
                                                               const std::vector<unsigned int> &n_phase_transitions_per_composition) const
      {
        // If strain rate is zero (like during the first time step) set it to some very small number
        // to prevent a division-by-zero, and a floating point exception.
        // Otherwise, calculate the square-root of the norm of the second invariant of the deviatoric-
        // strain rate (often simplified as epsilondot_ii)
        const double edot_ii = std::max(std::sqrt(std::max(-second_invariant(deviator(strain_rate)), 0.)),
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
        // minimum stress.
        const double damper_stress = 2. * damper_viscosity * (edot_ii - partial_strain_rates[4]);
        const double total_stress = viscoplastic_stress + damper_stress;

        // 6) Return the effective creep viscosity using the total stress
        return total_stress / (2. * edot_ii);
      }



      template <int dim>
      std::pair<double, double>
      CompositeViscoPlastic<dim>::calculate_isostress_log_strain_rate_and_derivative(const std::vector<std::array<std::pair<double, double>, 4>> &logarithmic_strain_rates_and_stress_derivatives,
                                                                                     const double viscoplastic_stress,
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
        const double f = viscoplastic_stress / (2. * maximum_viscosity);
        const double strain_rate = (strain_rate_scaling_factor * viscoplastic_strain_rate_sum) + f;
        partial_strain_rates[4] = strain_rate - viscoplastic_strain_rate_sum;
        // And the partial derivative of the log *total* strain rate
        // with respect to log *viscoplastic* stress follows as
        const double log_strain_rate_derivative = (strain_rate_scaling_factor * viscoplastic_strain_rate_sum * log_viscoplastic_strain_rate_derivative + f) / strain_rate;

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
                std::vector<double> partial_strain_rates_composition(n_decomposed_strain_rates, 0.);
                viscosity += (volume_fractions[composition]
                              * compute_composition_viscosity (pressure,
                                                               temperature,
                                                               grain_size,
                                                               composition,
                                                               strain_rate,
                                                               partial_strain_rates_composition,
                                                               phase_function_values,
                                                               n_phase_transitions_per_composition));
                for (unsigned int j=0; j < n_decomposed_strain_rates; ++j)
                  partial_strain_rates[j] += volume_fractions[composition] * partial_strain_rates_composition[j];
              }
            else
              {
                total_volume_fraction -= volume_fractions[composition];
              }
          }

        viscosity /= total_volume_fraction;

        for (unsigned int j=0; j < n_decomposed_strain_rates; ++j)
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
        // minimum stress.
        const double damper_stress = 2. * damper_viscosity * (edot_ii - partial_strain_rates[4]);
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
        partial_strain_rates[4] = strain_rate - viscoplastic_strain_rate_sum;
        // And the partial derivative of the log *total* strain rate
        // with respect to log *viscoplastic* stress follows as
        const double log_strain_rate_derivative = (strain_rate_scaling_factor * viscoplastic_strain_rate_sum * log_viscoplastic_strain_rate_derivative + f) / strain_rate;

        return std::make_pair(std::log(strain_rate), log_strain_rate_derivative);
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

        // Diffusion creep parameters
        Rheology::DiffusionCreep<dim>::declare_parameters(prm);

        // Dislocation creep parameters
        Rheology::DislocationCreep<dim>::declare_parameters(prm);

        // Dislocation creep parameters
        Rheology::PeierlsCreep<dim>::declare_parameters(prm);

        // Drucker Prager parameters
        Rheology::DruckerPragerPower<dim>::declare_parameters(prm);

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
        const double minimum_viscosity = prm.get_double("Minimum viscosity");

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
                    ExcMessage("You need to include at least one deformation mechanism."));

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
