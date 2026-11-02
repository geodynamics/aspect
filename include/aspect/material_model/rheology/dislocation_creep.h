// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2019-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_material_model_rheology_dislocation_creep_h
#define _aspect_material_model_rheology_dislocation_creep_h

#include <aspect/global.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      /**
       * Data structure for dislocation creep parameters.
       */
      struct DislocationCreepParameters
      {
        /**
         * The dislocation creep prefactor, activation energy, activation volume
         * and stress exponent.
         */
        double prefactor;
        double activation_energy;
        double activation_volume;
        double stress_exponent;

        /**
         * Constructor. Initializes all values to NaN.
         */
        DislocationCreepParameters();
      };

      template <int dim>
      class DislocationCreep : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          DislocationCreep();

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
           * Compute the creep parameters for the dislocation creep law.
           * If @p n_phase_transitions_per_composition points to a vector of
           * unsigned integers this is considered the number of phase transitions
           * for each compositional field and viscosity will be first computed on
           * each phase and then averaged for each compositional field.
           */
          const DislocationCreepParameters
          compute_creep_parameters (const unsigned int composition,
                                    const std::vector<double> &phase_function_values = {},
                                    const std::vector<unsigned int> &n_phase_transitions_per_composition = {}) const;

          /**
           * Compute the viscosity based on the dislocation creep law.
           * If @p n_phase_transitions_per_composition points to a vector of
           * unsigned integers this is considered the number of phase transitions
           * for each compositional field and viscosity will be first computed on
           * each phase and then averaged for each compositional field.
           */
          double
          compute_viscosity (const double strain_rate,
                             const double pressure,
                             const double temperature,
                             const unsigned int composition,
                             const std::vector<double> &phase_function_values = {},
                             const std::vector<unsigned int> &n_phase_transitions_per_composition = {}) const;

          /**
           * Compute the strain rate and first stress derivative
           * as a function of stress based on the dislocation creep law.
           */
          std::pair<double, double>
          compute_strain_rate_and_derivative (const double stress,
                                              const double pressure,
                                              const double temperature,
                                              const DislocationCreepParameters creep_parameters) const;

          /**
           * Compute the logarithm of strain rate and first derivative with respect to
           * the logarithm of the stress based on the dislocation creep law.
           */
          std::pair<double, double>
          compute_log_strain_rate_and_derivative (const double log_stress,
                                                  const double pressure,
                                                  const double temperature,
                                                  const DislocationCreepParameters creep_parameters) const;

        private:

          /**
           * List of dislocation creep prefactors A.
           *
           * This variable is read from the parameter file through a parameter called 'Prefactors for dislocation creep'.
           */
          std::vector<double> prefactors;

          /**
           * List of dislocation creep stress exponents n.
           * This variable is read from the parameter file through a parameter called 'Stress exponents for dislocation creep'.
           */
          std::vector<double> stress_exponents;

          /**
           * List of dislocation creep activation energies E.
           * This variable is read from the parameter file through a parameter called 'Activation energies for dislocation creep'.
           */
          std::vector<double> activation_energies;

          /**
           * List of dislocation creep activation volumes V.
           *
           * This variable is read from the parameter file through a parameter called 'Activation volumes for dislocation creep'.
           */
          std::vector<double> activation_volumes;

      };
    }
  }
}
#endif
