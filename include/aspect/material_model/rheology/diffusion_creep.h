/*
  Copyright (C) 2019 - 2024 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_rheology_diffusion_creep_h
#define _aspect_material_model_rheology_diffusion_creep_h

#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <aspect/material_model/utilities.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      /**
       * Data structure for diffusion creep parameters.
       */
      struct DiffusionCreepParameters
      {
        /**
         * The diffusion creep prefactor, activation energy, activation volume
         * and grain size exponent.
         */
        double prefactor;
        double activation_energy;
        double activation_volume;
        double stress_exponent;
        double grain_size_exponent;

        /**
         * Constructor. Initializes all values to NaN.
         */
        DiffusionCreepParameters();
      };

      template <int dim>
      class DiffusionCreep : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          DiffusionCreep();

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
           * Compute the creep parameters for the diffusion creep law.
           * If @p n_phase_transitions_per_composition points to a vector of
           * unsigned integers this is considered the number of phase transitions
           * for each compositional field and viscosity will be first computed on
           * each phase and then averaged for each compositional field.
           */
          const DiffusionCreepParameters
          compute_creep_parameters (const unsigned int composition,
                                    const std::vector<double> &phase_function_values = std::vector<double>(),
                                    const std::vector<unsigned int> &n_phase_transitions_per_composition = std::vector<unsigned int>()) const;

          /**
           * Compute the viscosity based on the diffusion creep law with
           * the fixed grain size given in the input file.
           * If @p n_phase_transitions_per_composition points to a vector of
           * unsigned integers this is considered the number of phase transitions
           * for each compositional field and viscosity will be first computed on
           * each phase and then averaged for each compositional field.
           */
          double
          compute_viscosity (const double pressure,
                             const double temperature,
                             const unsigned int composition,
                             const std::vector<double> &phase_function_values = std::vector<double>(),
                             const std::vector<unsigned int> &n_phase_transitions_per_composition = std::vector<unsigned int>()) const;

          /**
           * Compute the viscosity based on the diffusion creep law for the given @p grain_size.
           * If @p n_phase_transitions_per_composition points to a vector of
           * unsigned integers this is considered the number of phase transitions
           * for each compositional field and viscosity will be first computed on
           * each phase and then averaged for each compositional field.
           */
          double
          compute_viscosity (const double pressure,
                             const double temperature,
                             const double grain_size,
                             const unsigned int composition,
                             const std::vector<double> &phase_function_values = std::vector<double>(),
                             const std::vector<unsigned int> &n_phase_transitions_per_composition = std::vector<unsigned int>()) const;

          /**
           * Compute the strain rate and first stress derivative as a function
           * of stress based on the diffusion creep law with
           * the fixed grain size given in the input file.
           */
          std::pair<double, double>
          compute_strain_rate_and_derivative (const double stress,
                                              const double pressure,
                                              const double temperature,
                                              const DiffusionCreepParameters creep_parameters) const;

          /**
           * Compute the strain rate and first stress derivative as a function
           * of stress based on the diffusion creep law for the given @p grain_size.
           */
          std::pair<double, double>
          compute_strain_rate_and_derivative (const double stress,
                                              const double pressure,
                                              const double temperature,
                                              const double grain_size,
                                              const DiffusionCreepParameters creep_parameters) const;

          /**
           * Compute the logarithm of strain rate and first derivative with respect to
           * the logarithm of the stress based on the diffusion creep law
           * with the fixed grain size given in the input file.
           */
          std::pair<double, double>
          compute_log_strain_rate_and_derivative (const double log_stress,
                                                  const double pressure,
                                                  const double temperature,
                                                  const DiffusionCreepParameters creep_parameters) const;

          /**
           * Compute the logarithm of strain rate and first derivative with respect to
           * the logarithm of the stress based on the diffusion creep law
           * for the given @p grain_size.
           */
          std::pair<double, double>
          compute_log_strain_rate_and_derivative (const double log_stress,
                                                  const double pressure,
                                                  const double temperature,
                                                  const double grain_size,
                                                  const DiffusionCreepParameters creep_parameters) const;

        private:

          /**
           * List of diffusion creep prefactors A.
           */
          std::vector<double> prefactors_diffusion;

          /**
           * List of diffusion creep stress exponents n (usually = 1).
           */
          std::vector<double> stress_exponents_diffusion;

          /**
           * List of diffusion creep grain size exponents m.
           */
          std::vector<double> grain_size_exponents_diffusion;

          /**
           * List of diffusion creep activation energies E.
           */
          std::vector<double> activation_energies_diffusion;

          /**
           * List of diffusion creep activation volumes V.
           */
          std::vector<double> activation_volumes_diffusion;

          /**
           * Diffusion creep grain size d.  This is read from the
           * input file, and is only used by the functions that do
           * not take the grain size as additional argument.
           */
          double fixed_grain_size;
      };
    }
  }
}
#endif
