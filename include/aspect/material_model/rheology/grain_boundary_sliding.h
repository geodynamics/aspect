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

#ifndef _aspect_material_model_rheology_grain_boundary_sliding_h
#define _aspect_material_model_rheology_grain_boundary_sliding_h

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
       * Data structure for grain boundary sliding parameters.
       */
      struct GrainBoundarySlidingParameters
      {
        /**
         * The grain boundary sliding prefactor, activation energy, activation volume
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
        GrainBoundarySlidingParameters();
      };

      template <int dim>
      class GrainBoundarySliding : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          GrainBoundarySliding();

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
           * Compute the creep parameters for the grain boundary sliding law.
           * If @p n_phase_transitions_per_composition points to a vector of
           * unsigned integers this is considered the number of phase transitions
           * for each compositional field and viscosity will be first computed on
           * each phase and then averaged for each compositional field.
           */
          const GrainBoundarySlidingParameters
          compute_slide_parameters (const unsigned int composition,
                                    const std::vector<double> &phase_function_values = {},
                                    const std::vector<unsigned int> &n_phase_transitions_per_composition = {}) const;

          /**
           * Compute the viscosity based on the grain boundary sliding law with
           * the fixed grain size given in the input file.
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
           * Compute the viscosity based on the grain boundary sliding law for the given @p grain_size.
           * If @p n_phase_transitions_per_composition points to a vector of
           * unsigned integers this is considered the number of phase transitions
           * for each compositional field and viscosity will be first computed on
           * each phase and then averaged for each compositional field.
           */
          double
          compute_viscosity (const double strain_rate,
                             const double pressure,
                             const double temperature,
                             const double grain_size,
                             const unsigned int composition,
                             const std::vector<double> &phase_function_values = {},
                             const std::vector<unsigned int> &n_phase_transitions_per_composition = {}) const;

        private:

          /**
           * List of grain boundary sliding prefactors A.
           *
           * This variable is read from the parameter file through a parameter called 'Prefactors for grain boundary sliding'.
           */
          std::vector<double> prefactors;

          /**
           * List of grain boundary sliding stress exponents n (for ice = 1.8).
           * This variable is read from the parameter file through a parameter called 'Stress exponents for grain boundary sliding'.
           */
          std::vector<double> stress_exponents;

          /**
           * List of grain boundary sliding grain size exponents m.
           * This variable is read from the parameter file through a parameter called 'Grain size exponents for grain boundary sliding'.
           */
          std::vector<double> grain_size_exponents;

          /**
           * List of grain boundary sliding activation energies E.
           * This variable is read from the parameter file through a parameter called 'Activation energies for grain boundary sliding'.
           */
          std::vector<double> activation_energies;

          /**
           * List of grain boundary sliding activation volumes V.
           * This variable is read from the parameter file through a parameter called 'Activation volumes for grain boundary sliding'.
           */
          std::vector<double> activation_volumes;

          /**
           * Grain boundary sliding grain size d.  This is read from the
           * input file, and is only used by the functions that do
           * not take the grain size as additional argument.
           *
           * This variable is read from the parameter file through a parameter called 'Grain size'.
           */
          double fixed_grain_size;
      };
    }
  }
}
#endif
