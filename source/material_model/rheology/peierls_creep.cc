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


#include <aspect/material_model/rheology/peierls_creep.h>
#include <aspect/utilities.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      PeierlsCreepParameters::PeierlsCreepParameters()
        : prefactor (numbers::signaling_nan<double>()),
          stress_exponent (numbers::signaling_nan<double>()),
          activation_energy (numbers::signaling_nan<double>()),
          activation_volume (numbers::signaling_nan<double>()),
          peierls_stress (numbers::signaling_nan<double>()),
          glide_parameter_p (numbers::signaling_nan<double>()),
          glide_parameter_q (numbers::signaling_nan<double>()),
          fitting_parameter (numbers::signaling_nan<double>()),
          stress_cutoff (numbers::signaling_nan<double>())
      {}



      template <int dim>
      PeierlsCreep<dim>::PeierlsCreep ()
        = default;



      /**
       * Compute the creep parameters for the Peierls creep law.
       */
      template <int dim>
      const PeierlsCreepParameters
      PeierlsCreep<dim>::compute_creep_parameters (const unsigned int composition,
                                                   const std::vector<double> &phase_function_values,
                                                   const std::vector<unsigned int> &n_phase_transitions_per_composition) const
      {
        PeierlsCreepParameters creep_parameters;
        if (phase_function_values == std::vector<double>())
          {
            // no phases
            creep_parameters.prefactor = prefactors[composition];
            creep_parameters.stress_exponent = stress_exponents[composition];
            creep_parameters.activation_energy = activation_energies[composition];
            creep_parameters.activation_volume = activation_volumes[composition];
            creep_parameters.peierls_stress = peierls_stresses[composition];
            creep_parameters.glide_parameter_p = glide_parameters_p[composition];
            creep_parameters.glide_parameter_q = glide_parameters_q[composition];
            creep_parameters.fitting_parameter = fitting_parameters[composition];
            creep_parameters.stress_cutoff = stress_cutoffs[composition];
          }
        else
          {
            // Average among phases. This averaging is not strictly correct, but
            // it will not matter much if the parameters are similar across transitions.
            creep_parameters.prefactor = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                         prefactors, composition,  MaterialModel::MaterialUtilities::PhaseUtilities::logarithmic);
            creep_parameters.stress_exponent = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                               stress_exponents, composition);
            creep_parameters.activation_energy = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                                 activation_energies, composition);
            creep_parameters.activation_volume = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                                 activation_volumes, composition);
            creep_parameters.peierls_stress = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                              peierls_stresses, composition);
            creep_parameters.glide_parameter_p = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                                 glide_parameters_p, composition);
            creep_parameters.glide_parameter_q = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                                 glide_parameters_q, composition);
            creep_parameters.fitting_parameter = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                                 fitting_parameters, composition);
            creep_parameters.stress_cutoff = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                             stress_cutoffs, composition);
          }
        return creep_parameters;
      }



      template <int dim>
      double
      PeierlsCreep<dim>::compute_approximate_viscosity (const double strain_rate,
                                                        const double pressure,
                                                        const double temperature,
                                                        const unsigned int composition,
                                                        const std::vector<double> &phase_function_values,
                                                        const std::vector<unsigned int> &n_phase_transitions_per_composition) const
      {
        /**
         * An approximation of the Peierls creep formulation, where stress is replaced with strain rate
         * (second invariant) when calculating viscosity. This substitution requires a fitting parameter (gamma),
         * which describes the relationship between stress and the Peierls stress (gamma = stress/stress_peierls).
         * Details of this derivation can be found in the following document:
         * https://ucdavis.app.box.com/s/cl5mwhkjeabol4otrdukfcdwfvg9me4w/file/705438695737
         * The formulation for the viscosity is:
         *   stress_term * arrhenius_term * strain_rate_term, where
         * stress_term      = (0.5 * gamma * sigma_p) / (A * ((gamma * sigma_p)^n)^(1/(s+n)))
         * arrhenius_term   = exp(( (E + P * V) / (R * T)) * (1 - gamma^p)^q / (s + n) )
         * strain_rate_term = edot_ii^(1 / (s + n) - 1)
         * Above,
         * s = ((E + P * V) / (R * T)) * p * q * (1 - gamma^p)^(q-1) * (gamma^p)
         * sigma_p is the Peierls stress
         * gamma is the Peierls creep fitting parameter
         * p is the first Peierls glide parameter
         * q is the second Peierls creep glide parameter
         * A is the Peierls prefactor term (1/s 1/Pa^2)
         * n is the Peierls stress exponent
         * E is the Peierls activation energy (J/mol)
         * V is the Peierls activation volume (m^3/mol). However, to date V has always been set to 0.
         * edot_ii is the second invariant of the deviatoric strain rate (1/s)
         * T is temperature (K)
         * R is the gas constant
         */

        const PeierlsCreepParameters p = compute_creep_parameters(composition, phase_function_values, n_phase_transitions_per_composition);

        const double s = ( (p.activation_energy + pressure * p.activation_volume) / (constants::gas_constant * temperature)) *
                         p.glide_parameter_p * p.glide_parameter_q *
                         std::pow((1. - std::pow(p.fitting_parameter, p.glide_parameter_p)),(p.glide_parameter_q - 1.)) *
                         std::pow(p.fitting_parameter, p.glide_parameter_p);

        const double stress_term = 0.5 * (p.fitting_parameter * p.peierls_stress) /
                                   std::pow((p.prefactor * std::pow(p.fitting_parameter * p.peierls_stress,p.stress_exponent)),( 1. / (s + p.stress_exponent)));

        const double arrhenius_term = std::exp( ((p.activation_energy + pressure * p.activation_volume) / (constants::gas_constant * temperature)) *
                                                (std::pow((1. - std::pow(p.fitting_parameter,p.glide_parameter_p)),p.glide_parameter_q)) /
                                                (s + p.stress_exponent) );

        const double strain_rate_term = std::pow(strain_rate, (1. / (s + p.stress_exponent)) - 1.);

        double viscosity_peierls = stress_term * arrhenius_term * strain_rate_term;

        Assert (viscosity_peierls > 0.0,
                ExcMessage ("Negative peierls viscosity detected. This is unphysical and should not happen. "
                            "Check for negative parameters."));

        // Creep viscosities become extremely large at low
        // temperatures and can therefore provoke floating-point overflow errors. In
        // real rocks, other deformation mechanisms become dominant at low temperatures,
        // so these high viscosities are never achieved. It is therefore both reasonable
        // and desirable to require the single-mechanism viscosity to be smaller than
        // std::sqrt(max_double).
        viscosity_peierls = std::min(viscosity_peierls, std::sqrt(std::numeric_limits<double>::max()));

        return viscosity_peierls;
      }



      template <int dim>
      double
      PeierlsCreep<dim>::compute_exact_viscosity (const double strain_rate,
                                                  const double pressure,
                                                  const double temperature,
                                                  const unsigned int composition,
                                                  const std::vector<double> &phase_function_values,
                                                  const std::vector<unsigned int> &n_phase_transitions_per_composition) const
      {
        /**
         * A generalized Peierls creep formulation. The Peierls creep expression
         * for the strain rate has multiple stress-dependent terms, and cannot be
         * directly inverted to find an expression for viscosity in terms of
         * strain rate. For this reason, a Newton-Raphson iteration is required,
         * which can be quite expensive.
         * The equation for the strain rate is given in
         * compute_exact_strain_rate_and_derivative.
         */
        const PeierlsCreepParameters p = compute_creep_parameters(composition, phase_function_values, n_phase_transitions_per_composition);

        // The generalized Peierls creep flow law cannot be expressed as viscosity in
        // terms of strain rate, because there are two stress-dependent terms
        // in the strain rate expression.
        // We use Newton's method to find the second invariant of the stress tensor.

        // Apply a strict cutoff if this option is chosen by user. A strain rate cutoff
        // will be first computed and then compared to the input strain rate. A cutoff
        // on stress will be triggered if the input strain rate is smaller.
        const double log_strain_rate = std::log(strain_rate);

        if (apply_strict_cutoff)
          {
            const std::pair<double, double> log_edot_and_deriv = compute_exact_log_strain_rate_and_derivative(std::log(p.stress_cutoff), pressure, temperature, p);
            if (log_strain_rate < log_edot_and_deriv.first)
              {
                double viscosity = 0.5 * p.stress_cutoff / strain_rate;
                return viscosity;
              }
          }

        // Create a starting guess for the stress using
        // the approximate form of the viscosity expression
        double viscosity = compute_approximate_viscosity(strain_rate, pressure, temperature, composition);
        double log_stress_ii = std::log(2.*viscosity*strain_rate);

        // Before the first iteration, compute the residual
        // of the initial guess and the derivative
        unsigned int stress_iteration = 0;
        const std::pair<double, double> log_edot_and_deriv = compute_exact_log_strain_rate_and_derivative(log_stress_ii, pressure, temperature, p);
        double strain_rate_residual = log_edot_and_deriv.first - log_strain_rate;
        double log_strain_rate_deriv = log_edot_and_deriv.second;

        while (std::abs(strain_rate_residual) > strain_rate_residual_threshold
               && stress_iteration < stress_max_iteration_number)
          {
            // If the strain rate derivative is zero, we catch it below.
            if (log_strain_rate_deriv>std::numeric_limits<double>::min())
              log_stress_ii -= strain_rate_residual/log_strain_rate_deriv;

            const std::pair<double, double> log_edot_and_deriv = compute_exact_log_strain_rate_and_derivative(log_stress_ii, pressure, temperature, p);

            strain_rate_residual = log_edot_and_deriv.first - log_strain_rate;
            log_strain_rate_deriv = log_edot_and_deriv.second;

            stress_iteration += 1;

            // If anything that would be used in the next iteration is not finite, the
            // Newton iteration would trigger an exception and we want to abort the
            // iteration instead.
            // Currently, we still throw an exception, but if this exception is thrown,
            // another more robust iterative scheme should be implemented
            // (similar to that seen in the diffusion-dislocation material model).
            const bool abort_newton_iteration = !numbers::is_finite(log_stress_ii)
                                                || !numbers::is_finite(strain_rate_residual)
                                                || !numbers::is_finite(log_strain_rate_deriv)
                                                || log_strain_rate_deriv < std::numeric_limits<double>::min()
                                                || stress_iteration == stress_max_iteration_number;
            AssertThrow(!abort_newton_iteration,
                        ExcMessage("No convergence has been reached in the loop that determines "
                                   "the Peierls creep viscosity. Aborting! "
                                   "Residual is " + Utilities::to_string(strain_rate_residual) +
                                   " after " + Utilities::to_string(stress_iteration) + " iterations. "
                                   "You can increase the number of iterations by adapting the "
                                   "parameter 'Maximum Peierls strain rate iterations'."));
          }

        viscosity = 0.5*std::exp(log_stress_ii)/strain_rate;

        return viscosity;
      }



      template <int dim>
      double
      PeierlsCreep<dim>::compute_viscosity (const double strain_rate,
                                            const double pressure,
                                            const double temperature,
                                            const unsigned int composition,
                                            const std::vector<double> &phase_function_values,
                                            const std::vector<unsigned int> &n_phase_transitions_per_composition) const
      {
        double viscosity = 0.0;

        switch (peierls_creep_flow_law)
          {
            case viscosity_approximation:
            {
              viscosity = compute_approximate_viscosity(strain_rate, pressure, temperature, composition, phase_function_values, n_phase_transitions_per_composition);
              break;
            }
            case exact:
            {
              viscosity = compute_exact_viscosity(strain_rate, pressure, temperature, composition, phase_function_values, n_phase_transitions_per_composition);
              break;
            }
            default:
            {
              AssertThrow(false, ExcNotImplemented());
              break;
            }
          }

        return viscosity;
      }



      template <int dim>
      std::pair<double, double>
      PeierlsCreep<dim>::compute_approximate_strain_rate_and_derivative (const double stress,
                                                                         const double pressure,
                                                                         const double temperature,
                                                                         const PeierlsCreepParameters creep_parameters) const
      {
        /**
        * b = (E+P*V)/(R*T)
        * c = std::pow(gamma, p)
        * d = std::pow(1. - c, q)
        * s = b*p*q*c*d/(1. - c)
        * arrhenius = std::exp(-b*d)
        *
        * edot_ii = A * std::pow(stress, s + n) * std::pow(gamma*peierls_stress, -s) * arrhenius
        * deriv = edot_ii / stress * (s + n)
        */
        const PeierlsCreepParameters p = creep_parameters;

        const double b = (p.activation_energy + pressure*p.activation_volume)/(constants::gas_constant * temperature);
        const double c = std::pow(p.fitting_parameter, p.glide_parameter_p);
        const double d = std::pow(1. - c, p.glide_parameter_q);
        const double s = b*p.glide_parameter_p*p.glide_parameter_q*c*d/(1. - c);
        const double arrhenius = std::exp(-b*d);

        const double edot_ii = p.prefactor * std::pow(stress, s + p.stress_exponent) * std::pow(p.fitting_parameter*p.peierls_stress, -s) * arrhenius;
        const double deriv = edot_ii / stress * (s + p.stress_exponent);

        return std::make_pair(edot_ii, deriv);
      }



      template <int dim>
      std::pair<double, double>
      PeierlsCreep<dim>::compute_approximate_log_strain_rate_and_derivative (const double log_stress,
                                                                             const double pressure,
                                                                             const double temperature,
                                                                             const PeierlsCreepParameters creep_parameters) const
      {
        /**
        * b = (E+P*V)/(R*T)
        * c = std::pow(gamma, p)
        * d = std::pow(1. - c, q)
        * s = b*p*q*c*d/(1. - c)
        * log_arrhenius = -b*d
        *
        * log_strain_rate = log_A + (s + n) * std::log(stress) - s * std::log(gamma*peierls_stress) + log_arrhenius
        * total_stress_exponent = (s + n)
        */
        const PeierlsCreepParameters p = creep_parameters;

        const double b = (p.activation_energy + pressure*p.activation_volume)/(constants::gas_constant * temperature);
        const double c = std::pow(p.fitting_parameter, p.glide_parameter_p);
        const double d = std::pow(1. - c, p.glide_parameter_q);
        const double s = b*p.glide_parameter_p*p.glide_parameter_q*c*d/(1. - c);
        const double log_arrhenius = -b*d;

        const double total_stress_exponent = s + p.stress_exponent;
        const double log_strain_rate = std::log(p.prefactor) + total_stress_exponent * log_stress - s * std::log(p.fitting_parameter*p.peierls_stress) + log_arrhenius;

        return std::make_pair(log_strain_rate, total_stress_exponent);
      }



      template <int dim>
      std::pair<double, double>
      PeierlsCreep<dim>::compute_exact_strain_rate_and_derivative (const double stress,
                                                                   const double pressure,
                                                                   const double temperature,
                                                                   const PeierlsCreepParameters creep_parameters) const
      {
        /**
        * b = (E+P*V)/(R*T)
        * c = std::pow(stress/peierls_stress, p)
        * d = std::pow(1 - c, q)
        * s = b*p*q*c*d/(1 - c)
        * arrhenius = std::exp(-b*d)
        *
        * edot_ii = A * std::pow(stress, n) * arrhenius
        * deriv = edot_ii / stress * (s + n)
        */
        const PeierlsCreepParameters p = creep_parameters;
        if (stress < p.stress_cutoff)
          {

            /**
            * For Peierls creep flow laws that have a stress exponent equal to zero the strain rate does not approach zero as
            * stress approaches zero. To ensure convergence in the solver, the strain rate is modelled as a quadratic function
            * of stress;
            * edot_ii = quadratic_term*stress^2 + linear_term*stress
            * Where the quadratic and linear terms are defined at a constant cutoff temperature and pressure.
            * T_cutoff = (E/R), P_cutoff = 0
            * s_cutoff = p*q*c_cutoff*d_cutoff / (1 - c_cutoff)
            * arrhenius_cutoff = std::exp(-d_cutoff)
            */
            const double c_cutoff = std::pow(p.stress_cutoff/p.peierls_stress, p.glide_parameter_p);
            const double d_cutoff = std::pow(1. - c_cutoff, p.glide_parameter_q);
            const double s_cutoff = p.glide_parameter_p*p.glide_parameter_q*c_cutoff*d_cutoff/(1. - c_cutoff);
            const double arrhenius_cutoff = std::exp(-d_cutoff);
            const double edot_ii_cutoff = p.prefactor * std::pow(p.stress_cutoff, p.stress_exponent) * arrhenius_cutoff;
            const double deriv_cutoff = edot_ii_cutoff / p.stress_cutoff * (s_cutoff + p.stress_exponent);
            const double quadratic_term = (deriv_cutoff - edot_ii_cutoff / p.stress_cutoff) / p.stress_cutoff / arrhenius_cutoff;
            const double linear_term = (2*(edot_ii_cutoff / p.stress_cutoff) - deriv_cutoff) / arrhenius_cutoff;

            const double b = (p.activation_energy + pressure*p.activation_volume)/(constants::gas_constant * temperature);
            const double arrhenius = std::exp(-b*d_cutoff);
            const double edot_ii = (quadratic_term*Utilities::fixed_power<2>(stress) + linear_term*stress) * arrhenius;
            const double deriv = (2*quadratic_term*stress + linear_term) * arrhenius;

            return std::make_pair(edot_ii, deriv);
          }

        else
          {
            const double b = (p.activation_energy + pressure*p.activation_volume)/(constants::gas_constant * temperature);
            const double c = std::pow(stress/p.peierls_stress, p.glide_parameter_p);
            const double d = std::pow(1. - c, p.glide_parameter_q);
            const double s = b*p.glide_parameter_p*p.glide_parameter_q*c*d/(1. - c);
            const double arrhenius = std::exp(-b*d);

            const double edot_ii = p.prefactor * std::pow(stress, p.stress_exponent) * arrhenius;
            const double deriv = edot_ii / stress * (s + p.stress_exponent);

            return std::make_pair(edot_ii, deriv);
          }
      }



      template <int dim>
      std::pair<double, double>
      PeierlsCreep<dim>::compute_exact_log_strain_rate_and_derivative (const double log_stress,
                                                                       const double pressure,
                                                                       const double temperature,
                                                                       const PeierlsCreepParameters creep_parameters) const
      {
        /**
        * b = (E+P*V)/(R*T)
        * c = std::pow(stress/peierls_stress, p)
        * d = std::pow(1 - c, q)
        *
        * log_edot_ii = std::log(A) + n * std::log(stress) - b*d
        * deriv_log = n + p * q * b * std::(1-c, p.q - 1)
        * The deriv_log is the derivative of log(edot_ii) to log(stress).
        */
        const PeierlsCreepParameters p = creep_parameters;
        const double stress = std::exp(log_stress);
        if (stress < p.stress_cutoff)
          {

            /**
            * For Peierls creep flow laws that have a stress exponent equal to zero the strain rate does not approach zero as
            * stress approaches zero. To ensure convergence in the solver, the strain rate is modelled as a quadratic function
            * of stress;
            * edot_ii = quadratic_term*stress^2 + linear_term*stress
            * Where the quadratic and linear terms are defined at a constant cutoff temperature and pressure.
            * T_cutoff = (E/R), P_cutoff = 0
            * s_cutoff = p*q*c_cutoff*d_cutoff / (1 - c_cutoff)
            * arrhenius_cutoff = std::exp(-d_cutoff)
            */
            const double c_cutoff = std::pow(p.stress_cutoff/p.peierls_stress, p.glide_parameter_p);
            const double d_cutoff = std::pow(1. - c_cutoff, p.glide_parameter_q);
            const double s_cutoff = p.glide_parameter_p*p.glide_parameter_q*c_cutoff*d_cutoff/(1. - c_cutoff);
            const double arrhenius_cutoff = std::exp(-d_cutoff);
            const double edot_ii_cutoff = p.prefactor * std::pow(p.stress_cutoff, p.stress_exponent) * arrhenius_cutoff;
            const double deriv_cutoff = edot_ii_cutoff / p.stress_cutoff * (s_cutoff + p.stress_exponent);
            const double quadratic_term = (deriv_cutoff - edot_ii_cutoff / p.stress_cutoff) / p.stress_cutoff / arrhenius_cutoff;
            const double linear_term = (2*(edot_ii_cutoff / p.stress_cutoff) - deriv_cutoff) / arrhenius_cutoff;

            const double b = (p.activation_energy + pressure*p.activation_volume)/(constants::gas_constant * temperature);
            const double arrhenius = std::exp(-b*d_cutoff);
            const double edot_ii = (quadratic_term*Utilities::fixed_power<2>(stress) + linear_term*stress) * arrhenius;
            const double deriv_log = 2 - linear_term / (quadratic_term * stress + linear_term);

            return std::make_pair(std::log(edot_ii), deriv_log);
          }
        else
          {
            const double b = (p.activation_energy + pressure*p.activation_volume)/(constants::gas_constant * temperature);
            const double c = std::pow(stress/p.peierls_stress, p.glide_parameter_p);
            const double d = std::pow(1. - c, p.glide_parameter_q);

            const double log_edot_ii = std::log(p.prefactor) + p.stress_exponent * log_stress - b*d ;
            const double deriv_log = p.stress_exponent + p.glide_parameter_p * p.glide_parameter_q * b * c * std::pow(1-c, p.glide_parameter_q - 1);

            return std::make_pair(log_edot_ii, deriv_log);
          }
      }



      template <int dim>
      std::pair<double, double>
      PeierlsCreep<dim>::compute_strain_rate_and_derivative (const double stress,
                                                             const double pressure,
                                                             const double temperature,
                                                             const PeierlsCreepParameters creep_parameters) const
      {
        std::pair<double, double> strain_rate_and_deriv;

        switch (peierls_creep_flow_law)
          {
            case viscosity_approximation:
            {
              strain_rate_and_deriv = compute_approximate_strain_rate_and_derivative(stress, pressure, temperature, creep_parameters);
              break;
            }
            case exact:
            {
              strain_rate_and_deriv = compute_exact_strain_rate_and_derivative(stress, pressure, temperature, creep_parameters);
              break;
            }
            default:
            {
              AssertThrow(false, ExcNotImplemented());
              break;
            }
          }
        return strain_rate_and_deriv;
      }



      template <int dim>
      void
      PeierlsCreep<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Peierls creep flow law", "viscosity approximation",
                           Patterns::Selection("viscosity approximation|exact"),
                           "Select what type of Peierls creep flow law to use. Currently, the "
                           "available options are 'exact', which uses a Newton-Raphson iterative method "
                           "to find the stress and then compute viscosity, and 'viscosity approximation', "
                           "in which viscosity is an explicit function of the strain rate invariant, "
                           "rather than stress. ");

        // Viscosity iteration parameters
        prm.declare_entry ("Peierls strain rate residual tolerance", "1e-10", Patterns::Double(0.),
                           "Tolerance for the iterative solve to find the correct Peierls creep strain rate. "
                           "The tolerance is expressed as the difference between the natural logarithm of the "
                           "input strain rate and the strain rate at the current iteration.");
        prm.declare_entry ("Maximum Peierls strain rate iterations", "40", Patterns::Integer(0),
                           "Maximum number of iterations to find the correct "
                           "Peierls strain rate.");

        // Rheological parameters
        prm.declare_entry ("Prefactors for Peierls creep", "1.4e-19",
                           Patterns::Anything(),
                           "List of viscosity prefactors, $A$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\pascal}$^{-n_{\\text{peierls}}}$ \\si{\\per\\second}");
        prm.declare_entry ("Stress exponents for Peierls creep", "2.0",
                           Patterns::Anything(),
                           "List of stress exponents, $n_{\\text{peierls}}$, for background material and compositional "
                           "fields, for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value.  Units: None.");
        prm.declare_entry ("Activation energies for Peierls creep", "320e3",
                           Patterns::Anything(),
                           "List of activation energies, $E$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value. Units: \\si{\\joule\\per\\mole}.");
        prm.declare_entry ("Activation volumes for Peierls creep", "1.4e-5",
                           Patterns::Anything(),
                           "List of activation volumes, $V$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\meter\\cubed\\per\\mole}.");
        prm.declare_entry ("Peierls stresses", "5.e9",
                           Patterns::Anything(),
                           "List of stress limits for Peierls creep $\\sigma_{\\text{peierls}}$ for background "
                           "material and compositional fields, for a total of N+1 values, where N is the number "
                           "of all compositional fields or only those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\pascal}");
        prm.declare_entry ("Peierls fitting parameters", "0.17",
                           Patterns::Anything(),
                           "List of fitting parameters $\\gamma$ between stress $\\sigma$ and the Peierls "
                           "stress $\\sigma_{\\text{peierls}}$ for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. If only one value is given, "
                           "then all use the same value. Units: none");
        prm.declare_entry ("Peierls glide parameters p", "0.5",
                           Patterns::Anything(),
                           "List of the first Peierls creep glide parameters, $p$, for background and compositional "
                           "fields for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value. Units: none");
        prm.declare_entry ("Peierls glide parameters q", "1.0",
                           Patterns::Anything(),
                           "List of the second Peierls creep glide parameters, $q$, for background and compositional "
                           "fields for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "If only one value is given, then all use the same value. Units: none");
        prm.declare_entry ("Cutoff stresses for Peierls creep", "0.0",
                           Patterns::Anything(),
                           "List of the Stress thresholds below which the strain rate is solved for as a quadratic "
                           "function of stress to aid with convergence when stress exponent n=0. Units: \\si{\\pascal}");
        prm.declare_entry ("Apply strict stress cutoff for Peierls creep", "false", Patterns::Bool(),
                           "Whether the cutoff stresses for Peierls creep are used as the minimum "
                           "stresses in the Peierls rheology");

      }



      template <int dim>
      void
      PeierlsCreep<dim>::parse_parameters (ParameterHandler &prm,
                                           const std::unique_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition)
      {
        // Retrieve the list of composition names
        std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();

        // Retrieve the list of names of fields that represent chemical compositions, and not, e.g.,
        // plastic strain
        std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();

        // Establish that a background field is required here
        compositional_field_names.insert(compositional_field_names.begin(), "background");
        chemical_field_names.insert(chemical_field_names.begin(),"background");

        if (prm.get ("Peierls creep flow law") == "viscosity approximation")
          peierls_creep_flow_law = viscosity_approximation;
        else if (prm.get ("Peierls creep flow law") == "exact")
          peierls_creep_flow_law = exact;
        else
          AssertThrow(false, ExcMessage("Not a valid Peierls creep flow law"));

        // Iteration parameters
        strain_rate_residual_threshold = prm.get_double ("Peierls strain rate residual tolerance");
        stress_max_iteration_number = prm.get_integer ("Maximum Peierls strain rate iterations");

        // Rheological parameters
        // Make options file for parsing maps to double arrays
        Utilities::MapParsing::Options options(chemical_field_names, "Prefactors for Peierls creep");
        options.list_of_allowed_keys = compositional_field_names;
        options.allow_multiple_values_per_key = true;
        if (expected_n_phases_per_composition)
          {
            options.n_values_per_key = *expected_n_phases_per_composition;

            // check_values_per_key is required to be true to duplicate single values
            // if they are to be used for all phases associated with a given key.
            options.check_values_per_key = true;
          }

        prefactors = Utilities::MapParsing::parse_map_to_double_array(prm.get("Prefactors for Peierls creep"),
                                                                      options);

        options.property_name = "Stress exponents for Peierls creep";
        stress_exponents = Utilities::MapParsing::parse_map_to_double_array(prm.get("Stress exponents for Peierls creep"),
                                                                            options);

        options.property_name = "Activation energies for Peierls creep";
        activation_energies = Utilities::MapParsing::parse_map_to_double_array(prm.get("Activation energies for Peierls creep"),
                                                                               options);

        options.property_name = "Activation volumes for Peierls creep";
        activation_volumes = Utilities::MapParsing::parse_map_to_double_array(prm.get("Activation volumes for Peierls creep"),
                                                                              options);

        options.property_name = "Peierls stresses";
        peierls_stresses = Utilities::MapParsing::parse_map_to_double_array(prm.get("Peierls stresses"),
                                                                            options);

        options.property_name = "Peierls fitting parameters";
        fitting_parameters = Utilities::MapParsing::parse_map_to_double_array(prm.get("Peierls fitting parameters"),
                                                                              options);

        options.property_name = "Peierls glide parameters p";
        glide_parameters_p = Utilities::MapParsing::parse_map_to_double_array(prm.get("Peierls glide parameters p"),
                                                                              options);

        options.property_name = "Peierls glide parameters q";
        glide_parameters_q = Utilities::MapParsing::parse_map_to_double_array(prm.get("Peierls glide parameters q"),
                                                                              options);

        options.property_name = "Cutoff stresses for Peierls creep";
        stress_cutoffs = Utilities::MapParsing::parse_map_to_double_array(prm.get("Cutoff stresses for Peierls creep"),
                                                                          options);

        apply_strict_cutoff = prm.get_bool("Apply strict stress cutoff for Peierls creep");
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
    template class PeierlsCreep<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
