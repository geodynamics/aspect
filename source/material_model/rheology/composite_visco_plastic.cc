/*
  Copyright (C) 2020 - 2023 by the authors of the ASPECT code.

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


      template <int dim>
      CompositeViscoPlastic<dim>::CompositeViscoPlastic ()
        = default;



      template <int dim>
      double
      CompositeViscoPlastic<dim>::compute_viscosity (const double pressure,
                                                     const double temperature,
                                                     const std::vector<double> &volume_fractions,
                                                     const SymmetricTensor<2,dim> &strain_rate,
                                                     std::vector<double> &partial_strain_rates,
                                                     const std::vector<double> &phase_function_values,
                                                     const std::vector<unsigned int> &n_phase_transitions_per_composition) const
      {
        double viscosity = 0.;
        partial_strain_rates.resize(5, 0.);

        // Isostrain averaging
        double total_volume_fraction = 1.;
        for (unsigned int composition=0; composition < number_of_compositions; ++composition)
          {
            // Only include the contribution to the viscosity
            // from a given composition if the volume fraction exceeds
            // a certain (small) fraction.
            if (volume_fractions[composition] > 2.*std::numeric_limits<double>::epsilon())
              {
                std::vector<double> partial_strain_rates_composition(5, 0.);
                viscosity += (volume_fractions[composition]
                              * compute_composition_viscosity (pressure,
                                                               temperature,
                                                               composition,
                                                               strain_rate,
                                                               partial_strain_rates_composition,
                                                               phase_function_values,
                                                               n_phase_transitions_per_composition));
                for (unsigned int j=0; j < 5; ++j)
                  partial_strain_rates[j] += volume_fractions[composition] * partial_strain_rates_composition[j];
              }
            else
              {
                total_volume_fraction -= volume_fractions[composition];
              }

            viscosity /= total_volume_fraction;
            for (unsigned int j=0; j < 5; ++j)
              partial_strain_rates[j] /= total_volume_fraction;
          }
        return viscosity;
      }



      template <int dim>
      double
      CompositeViscoPlastic<dim>::compute_composition_viscosity (const double pressure,
                                                                 const double temperature,
                                                                 const unsigned int composition,
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
                                        min_strain_rate);
        const double log_edot_ii = std::log(edot_ii);

        Rheology::DiffusionCreepParameters diffusion_creep_parameters;
        Rheology::DislocationCreepParameters dislocation_creep_parameters;
        Rheology::PeierlsCreepParameters peierls_creep_parameters;
        Rheology::DruckerPragerParameters drucker_prager_parameters;

        // 1) Estimate the stress running through the creep elements and the maximum viscosity element.
        // These are all arranged in series. Taking the minimum viscosity assuming that all the strain
        // runs through a single element provides an excellent first approximation to the true viscosity.
        // The stress can then be calculated as 2 * eta * edot_ii
        double eta_creep_guess = limiting_creep_viscosity;

        if (use_diffusion_creep)
          {
            diffusion_creep_parameters = diffusion_creep->compute_creep_parameters(composition, phase_function_values, n_phase_transitions_per_composition);
            eta_creep_guess = std::min(diffusion_creep->compute_viscosity(pressure, temperature, composition, phase_function_values, n_phase_transitions_per_composition), eta_creep_guess);
          }

        if (use_dislocation_creep)
          {
            dislocation_creep_parameters = dislocation_creep->compute_creep_parameters(composition, phase_function_values, n_phase_transitions_per_composition);
            eta_creep_guess = std::min(dislocation_creep->compute_viscosity(edot_ii, pressure, temperature, composition, phase_function_values, n_phase_transitions_per_composition), eta_creep_guess);
          }

        if (use_peierls_creep)
          {
            peierls_creep_parameters = peierls_creep->compute_creep_parameters(composition);
            eta_creep_guess = std::min(peierls_creep->compute_approximate_viscosity(edot_ii, pressure, temperature, composition), eta_creep_guess);
          }

        if (use_drucker_prager)
          {
            drucker_prager_parameters = drucker_prager->compute_drucker_prager_parameters(composition, phase_function_values, n_phase_transitions_per_composition);
            eta_creep_guess = std::min(drucker_prager->compute_viscosity(drucker_prager_parameters.cohesion,
                                                                         drucker_prager_parameters.angle_internal_friction, pressure, edot_ii, drucker_prager_parameters.max_yield_stress), eta_creep_guess);
          }

        double log_creep_stress = std::log(2. * eta_creep_guess * edot_ii);

        // 2) Calculate the total strain rate given by this creep stress.
        // This will not be equal to the real total strain rate because
        // the creep stress has been calculated assuming that only one mechanism can
        // accommodate strain, whereas now we allow all the mechanisms to
        // accommodate strain.
        // Some strain is also partitioned into a "hard" viscosity damper,
        // which is designed to limit the effective viscosity between minimum and
        // maximum bounds and has an effective viscosity of eta_max - eta_min.
        std::array<std::pair<double, double>, 5> log_edot_and_deriv;

        if (use_diffusion_creep)
          log_edot_and_deriv[0] = diffusion_creep->compute_log_strain_rate_and_derivative(log_creep_stress, pressure, temperature, diffusion_creep_parameters);

        if (use_dislocation_creep)
          log_edot_and_deriv[1] = dislocation_creep->compute_log_strain_rate_and_derivative (log_creep_stress, pressure, temperature, dislocation_creep_parameters);

        if (use_peierls_creep)
          log_edot_and_deriv[2] = peierls_creep->compute_approximate_log_strain_rate_and_derivative(log_creep_stress, pressure, temperature, peierls_creep_parameters);

        if (use_drucker_prager)
          log_edot_and_deriv[3] = drucker_prager->compute_log_strain_rate_and_derivative (log_creep_stress, pressure, drucker_prager_parameters);

        log_edot_and_deriv[4] = std::make_pair(log_creep_stress - std::log(2. * limiting_creep_viscosity), 1.);

        // 3) Calculate the total log strain rate
        std::pair<double, double> log_edot_ii_and_deriv_iterate = calculate_log_strain_rate_and_deriv(log_edot_and_deriv);
        double strain_rate_residual = log_edot_ii_and_deriv_iterate.first - log_edot_ii;

        // In this rheology model, the total strain rate is partitioned between
        // different flow laws. We do not know how the strain is partitioned
        // between these flow laws, nor do we know the creep stress, which is
        // required to calculate the partitioning.

        // The following while loop conducts a Newton iteration to obtain the
        // creep stress, which we need in order to calculate the viscosity.
        unsigned int stress_iteration = 0;
        while (std::abs(strain_rate_residual) > strain_rate_residual_threshold
               && stress_iteration < stress_max_iteration_number)
          {
            // Apply the Newton update for the log creep stress using the
            // strain-rate residual and strain-rate stress derivative
            double delta_log_creep_stress = strain_rate_residual/log_edot_ii_and_deriv_iterate.second;
            log_creep_stress += delta_log_creep_stress;

            // Update the strain rates of all mechanisms with the new stress
            for (auto &i : active_flow_mechanisms)
              log_edot_and_deriv[i].first += log_edot_and_deriv[i].second * delta_log_creep_stress;

            // Compute the new log strain rate residual and log stress derivative
            log_edot_ii_and_deriv_iterate = calculate_log_strain_rate_and_deriv(log_edot_and_deriv);
            strain_rate_residual = log_edot_ii - log_edot_ii_and_deriv_iterate.first;

            ++stress_iteration;

            // If anything that would be used in the next iteration is not finite, the
            // Newton iteration would trigger an exception and we want to abort the
            // iteration instead.
            // Currently, we still throw an exception, but if this exception is thrown,
            // another more robust iterative scheme should be implemented
            // (similar to that seen in the diffusion-dislocation material model).
            const bool abort_newton_iteration = !numbers::is_finite(log_creep_stress)
                                                || !numbers::is_finite(strain_rate_residual)
                                                || stress_iteration == stress_max_iteration_number;
            AssertThrow(!abort_newton_iteration,
                        ExcMessage("No convergence has been reached in the loop that determines "
                                   "the composite viscous creep stress. Aborting! "
                                   "Residual is " + Utilities::to_string(strain_rate_residual) +
                                   " after " + Utilities::to_string(stress_iteration) + " iterations. "
                                   "You can increase the number of iterations by adapting the "
                                   "parameter 'Maximum creep strain rate iterations'."));
          }

        // The creep stress is not the total stress, so we still need to do a little work to obtain the effective viscosity.
        // First, we compute the stress running through the strain rate limiter, and then add that to the creep stress
        // NOTE: The viscosity of the strain rate limiter is equal to (minimum_viscosity*maximum_viscosity)/(maximum_viscosity - minimum_viscosity)
        const double creep_stress = std::exp(log_creep_stress);
        const double lim_stress = 2. * minimum_viscosity * (edot_ii - creep_stress / (2. * maximum_viscosity));
        const double total_stress = creep_stress + lim_stress;

        // Compute the strain rate experienced by the different mechanisms
        // These should sum to the total strain rate

        // The components of partial_strain_rates must be provided in the order
        // dictated by make_strain_rate_additional_outputs_names
        for (auto &i : active_flow_mechanisms)
          partial_strain_rates[i] = std::exp(log_edot_and_deriv[i].first);

        // Now we return the viscosity using the total stress
        return total_stress/(2.*edot_ii);
      }



      // Overload the + operator to act on two pairs of doubles.
      std::pair<double,double> operator+(const std::pair<double,double> &x, const std::pair<double,double> &y)
      {
        return std::make_pair(x.first+y.first, x.second+y.second);
      }



      template <int dim>
      std::pair<double, double>
      CompositeViscoPlastic<dim>::calculate_log_strain_rate_and_deriv(const std::array<std::pair<double, double>, 5> &logarithmic_strain_rates_and_stress_derivatives) const
      {
        // The total strain rate
        double strain_rate_sum = 0.0;

        // The sum of the stress derivatives multiplied by the mechanism strain rates
        double weighted_stress_derivative_sum = 0.0;

        // The first derivative of log(strain rate) with respect to log(stress)
        // is computed as sum_i(stress_exponent_i * edot_i) / sum_i(edot_i)
        // i.e., the stress exponents weighted by strain rate fraction
        // summed over the individual flow mechanisms (i).

        // Loop over all flow laws and add their contributions
        // to the strain rate and stress derivative
        for (const auto &log_strain_rate_and_stress_derivative : logarithmic_strain_rates_and_stress_derivatives)
          {
            double mechanism_log_strain_rate = log_strain_rate_and_stress_derivative.first;

            // Check if the exponent is within bounds to prevent underflow
            if (mechanism_log_strain_rate >= logmin)
              {
                const double mechanism_strain_rate = std::exp(mechanism_log_strain_rate);
                const double log_stress_derivative = log_strain_rate_and_stress_derivative.second;
                strain_rate_sum += mechanism_strain_rate;
                weighted_stress_derivative_sum += log_stress_derivative * mechanism_strain_rate;
              }
          }
        return std::make_pair(std::log(strain_rate_sum), weighted_stress_derivative_sum/strain_rate_sum);
      }



      template <int dim>
      void
      CompositeViscoPlastic<dim>::declare_parameters (ParameterHandler &prm)
      {
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
                           "creep strain rate.");

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
        // Retrieve the list of composition names
        const std::vector<std::string> list_of_composition_names = this->introspection().get_composition_names();

        // A background field is required by the subordinate material models
        number_of_compositions = list_of_composition_names.size() + 1;

        min_strain_rate = prm.get_double("Minimum strain rate");

        // Iteration parameters
        strain_rate_residual_threshold = prm.get_double ("Strain rate residual tolerance");
        stress_max_iteration_number = prm.get_integer ("Maximum creep strain rate iterations");

        // Read min and max viscosity parameters
        minimum_viscosity = prm.get_double ("Minimum viscosity");
        maximum_viscosity = prm.get_double ("Maximum viscosity");
        limiting_creep_viscosity = maximum_viscosity - minimum_viscosity;

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

        active_flow_mechanisms.push_back(4);

        AssertThrow(active_flow_mechanisms.size() != 1,
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
