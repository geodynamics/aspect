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

#ifndef _aspect_material_model_rheology_composite_visco_plastic_h
#define _aspect_material_model_rheology_composite_visco_plastic_h

#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <aspect/material_model/utilities.h>
#include <aspect/material_model/rheology/diffusion_creep.h>
#include <aspect/material_model/rheology/dislocation_creep.h>
#include <aspect/material_model/rheology/peierls_creep.h>
#include <aspect/material_model/rheology/drucker_prager_power.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      /**
       * A namespace for selecting how to average the viscosity
       * over multiple chemical compositions.
       *
       * Current options are:
       *  'isostress': Reuss (harmonic) averaging of viscosities.
       *  'isostrain': Voigt (arithmetic) averaging of viscosities.
       */
      namespace ViscosityAveraging
      {
        enum Kind
        {
          isostress,
          isostrain
        };

        /**
         * Read the viscosity averaging scheme from the parameter file
         * using the parameter name given in @p parameter_name and return
         * the enum that corresponds to this operation.
         */
        ViscosityAveraging::Kind
        parse (const std::string &parameter_name,
               const ParameterHandler &prm);
      }

      template <int dim>
      class CompositeViscoPlastic : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          CompositeViscoPlastic();

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
           * Compute the viscosity based on the composite viscous creep law.
           * If @p n_phase_transitions_per_composition points to a vector of
           * unsigned integers this is considered the number of phase transitions
           * for each compositional field and viscosity will be first computed on
           * each phase and then averaged for each compositional field.
           */
          double
          compute_viscosity (const double pressure,
                             const double temperature,
                             const double grain_size,
                             const std::vector<double> &volume_fractions,
                             const SymmetricTensor<2,dim> &strain_rate,
                             std::vector<double> &partial_strain_rates,
                             const std::vector<double> &phase_function_values = std::vector<double>(),
                             const std::vector<unsigned int> &n_phase_transitions_per_composition = std::vector<unsigned int>()) const;

        private:
          /**
           * Compute the isostress viscosity over all compositional fields
           * based on the composite viscous creep law.
           * If @p expected_n_phases_per_composition points to a vector of
           * unsigned integers this is considered the number of phase transitions
           * for each compositional field and viscosity will be first computed on
           * each phase and then averaged for each compositional field.
           */
          double
          compute_isostress_viscosity (const double pressure,
                                       const double temperature,
                                       const double grain_size,
                                       const std::vector<double> &volume_fractions,
                                       const SymmetricTensor<2,dim> &strain_rate,
                                       std::vector<double> &partial_strain_rates,
                                       const std::vector<double> &phase_function_values = std::vector<double>(),
                                       const std::vector<unsigned int> &n_phase_transitions_per_composition = std::vector<unsigned int>()) const;

          /**
           * Compute the total strain rate and the first derivative of log strain rate
           * with respect to log viscoplastic stress from individual log strain rate components
           * over all compositional fields. Also updates the partial_strain_rates vector.
           */
          std::pair<double, double>
          calculate_isostress_log_strain_rate_and_derivative(const std::vector<std::array<std::pair<double, double>, 4>> &logarithmic_strain_rates_and_stress_derivatives,
                                                             const double viscoplastic_stress,
                                                             std::vector<double> &partial_strain_rates) const;

          /**
           * Compute the isostrain viscosity over all compositional fields
           * based on the composite viscous creep law.
           * If @p expected_n_phases_per_composition points to a vector of
           * unsigned integers this is considered the number of phase transitions
           * for each compositional field and viscosity will be first computed on
           * each phase and then averaged for each compositional field.
           */
          double
          compute_isostrain_viscosity (const double pressure,
                                       const double temperature,
                                       const double grain_size,
                                       const std::vector<double> &volume_fractions,
                                       const SymmetricTensor<2,dim> &strain_rate,
                                       std::vector<double> &partial_strain_rates,
                                       const std::vector<double> &phase_function_values = std::vector<double>(),
                                       const std::vector<unsigned int> &n_phase_transitions_per_composition = std::vector<unsigned int>()) const;


          /**
           * Compute the viscosity for a single compositional field
           * based on the composite viscous creep law.
           * If @p n_phase_transitions_per_composition points to a vector of
           * unsigned integers this is considered the number of phase transitions
           * for each compositional field and viscosity will be first computed on
           * each phase and then averaged for each compositional field.
           */
          double
          compute_composition_viscosity (const double pressure,
                                         const double temperature,
                                         const double grain_size,
                                         const unsigned int composition,
                                         const SymmetricTensor<2,dim> &strain_rate,
                                         std::vector<double> &partial_strain_rates,
                                         const std::vector<double> &phase_function_values = std::vector<double>(),
                                         const std::vector<unsigned int> &n_phase_transitions_per_composition = std::vector<unsigned int>()) const;

          /**
           * Compute the total strain rate and the first derivative of log strain rate
           * with respect to log viscoplastic stress from individual log strain rate components
           * for a single compositional field.
           * Also updates the partial_strain_rates vector.
           */
          std::pair<double, double>
          calculate_composition_log_strain_rate_and_derivative(const std::array<std::pair<double, double>, 4> &logarithmic_strain_rates_and_stress_derivatives,
                                                               const double viscoplastic_stress,
                                                               std::vector<double> &partial_strain_rates) const;

          /**
           * Enumeration for selecting which type of viscosity averaging to use.
           */
          ViscosityAveraging::Kind viscosity_averaging_scheme;

          /**
           * Whether to use different deformation mechanisms.
           */
          bool use_diffusion_creep;
          bool use_dislocation_creep;
          bool use_peierls_creep;
          bool use_drucker_prager;

          /**
           * Vector storing which flow mechanisms are active.
           */
          std::vector<unsigned int> active_flow_mechanisms;

          /**
           * Total number of decomposed strain rates. This includes the maximum
           * viscosity limiter but not the minimum viscosity limiter,
           * which is arranged in parallel with the viscoplastic elements and
           * therefore does not contribute to the total strain rate.
           */
          static constexpr unsigned int n_decomposed_strain_rates = 5;

          /**
           * Pointers to objects for computing deformation mechanism
           * strain rates and effective viscosities.
           */
          std::unique_ptr<Rheology::DiffusionCreep<dim>> diffusion_creep;
          std::unique_ptr<Rheology::DislocationCreep<dim>> dislocation_creep;
          std::unique_ptr<Rheology::PeierlsCreep<dim>> peierls_creep;
          std::unique_ptr<Rheology::DruckerPragerPower<dim>> drucker_prager;

          /**
           * The expected number of chemical compositional fields.
           */
          unsigned int number_of_chemical_compositions;

          /**
           * The maximum viscosity, imposed via an isoviscous damper
           * in series with the composite viscoplastic element.
           */
          double maximum_viscosity;


          /**
           * The viscosity of an isoviscous damper placed in parallel
           * with the flow law components.
           * The viscosity is equal to the product of the
           * minimum and maximum viscosities
           * divided by the difference between the maximum and
           * minimum viscosities.
           */
          double damper_viscosity;

          /**
           * A factor used to scale the viscoplastic strain rate
           * up to the total strain rate.
           * The scaling factor is equal to the maximum viscosity
           * divided by the difference between the maximum and
           * minimum viscosities.
           */
          double strain_rate_scaling_factor;

          /**
           * The minimum strain rate allowed by the rheology.
           */
          double minimum_strain_rate;

          /**
           * The log strain rate threshold which must be passed
           * before successful termination of the Newton iteration
           * to determine the creep stress.
           */
          double log_strain_rate_residual_threshold;

          /**
           * The maximum number of iterations allowed to determine
           * the creep stress.
           */
          unsigned int stress_max_iteration_number;

          /**
           * Useful number to aid in adding together exponentials.
           */
          const double logmin = std::log(std::numeric_limits<double>::min());
      };
    }
  }
}
#endif
