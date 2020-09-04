/*
  Copyright (C) 2020 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_rheology_peierls_creep_h
#define _aspect_material_model_rheology_peierls_creep_h

#include <aspect/global.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    namespace Rheology
    {
      /**
       * Data structure for Peierls creep parameters.
       */
      struct PeierlsCreepParameters
      {
        /**
         * The Peierls creep prefactor, stress exponent, activation energy,
         * activation volume, Peierls stress, glide parameters p and q
         * and fitting parameter,
         */
        double prefactor;
        double stress_exponent;
        double activation_energy;
        double activation_volume;
        double peierls_stress;
        double glide_parameter_p;
        double glide_parameter_q;
        double fitting_parameter;
      };

      /**
       * Peierls creep is a low temperature, high stress viscous deformation mechanism.
       * Crystal deformation occurs through dislocation glide, which requires high stresses.
       * It is considered an important deformation mechanism in the lithosphere and subducting
       * slabs (Kumamoto et al., 2017, Science Advances). The approximate form of the Peierls
       * flow law used here is based on a derivation from Kameyama et al., 1999, Earth and
       * Planetary Science Letters.
       */
      template <int dim>
      class PeierlsCreep : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          PeierlsCreep();

          /**
           * Compute the creep parameters for the Peierls creep law.
           */
          const PeierlsCreepParameters
          compute_creep_parameters (const unsigned int composition,
                                    const std::vector<double> &phase_function_values = std::vector<double>(),
                                    const std::vector<unsigned int> &n_phases_per_composition = std::vector<unsigned int>()) const;

          /**
           * Declare the parameters this function takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters from the parameter file.
           * If @p expected_n_phases_per_composition points to a vector of
           * unsigned integers this is considered the number of phase transitions
           * for each compositional field and will be checked against the parsed
           * parameters.
           */
          void
          parse_parameters (ParameterHandler &prm,
                            const std::shared_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition =
                              std::shared_ptr<std::vector<unsigned int>>());

          /**
           * Compute the viscosity based on the approximate Peierls creep flow law.
           */
          double
          compute_approximate_viscosity (const double strain_rate,
                                         const double pressure,
                                         const double temperature,
                                         const unsigned int composition,
                                         const std::vector<double> &phase_function_values = std::vector<double>(),
                                         const std::vector<unsigned int> &n_phases_per_composition = std::vector<unsigned int>()) const;

          /**
           * Compute the viscosity based on the exact Peierls creep flow law.
           */
          double
          compute_exact_viscosity (const double strain_rate,
                                   const double pressure,
                                   const double temperature,
                                   const unsigned int composition,
                                   const std::vector<double> &phase_function_values = std::vector<double>(),
                                   const std::vector<unsigned int> &n_phases_per_composition = std::vector<unsigned int>()) const;

          /**
           * Compute the viscosity based on the selected Peierls creep flow law.
           * This function uses either the compute_approximate_viscosity
           * or the compute_exact_viscosity function.
           */
          double
          compute_viscosity (const double strain_rate,
                             const double pressure,
                             const double temperature,
                             const unsigned int composition,
                             const std::vector<double> &phase_function_values = std::vector<double>(),
                             const std::vector<unsigned int> &n_phases_per_composition = std::vector<unsigned int>()) const;

          /**
           * Compute the strain rate and first stress derivative
           * as a function of stress based on the approximate Peierls creep law.
           */
          std::pair<double, double>
          compute_approximate_strain_rate_and_derivative (const double stress,
                                                          const double pressure,
                                                          const double temperature,
                                                          const PeierlsCreepParameters creep_parameters) const;

          /**
           * Compute the strain rate and first stress derivative
           * as a function of stress based on the exact Peierls creep law.
           */
          std::pair<double, double>
          compute_exact_strain_rate_and_derivative (const double stress,
                                                    const double pressure,
                                                    const double temperature,
                                                    const PeierlsCreepParameters creep_parameters) const;

          /**
           * Compute the strain rate and first stress derivative
           * as a function of stress based on the selected Peierls creep law.
           * This function uses either the
           * compute_approximate_strain_rate_and_derivative
           * or the compute_exact_strain_rate_and_derivative function.
           */
          std::pair<double, double>
          compute_strain_rate_and_derivative (const double stress,
                                              const double pressure,
                                              const double temperature,
                                              const PeierlsCreepParameters creep_parameters) const;

        private:
          /**
           * Enumeration for selecting which type of Peierls creep flow law to use.
           * The options are currently "exact" (which requires a Newton-Raphson iteration
           * in compute_viscosity) and "viscosity_approximation", which
           * uses a simplified expression where the strain rate has a power law dependence
           * on the stress (so that the equation can be directly inverted for viscosity).
           * This approximation requires specifying one fitting parameter (gamma) to
           * obtain the best fit in the expected stress range for a given problem
           * (gamma = stress / peierls_stress). The derivation for this approximation
           * (derived by Magali Billen) can be found at:
           * https://ucdavis.app.box.com/s/cl5mwhkjeabol4otrdukfcdwfvg9me4w/file/705438695737.
           */
          enum PeierlsCreepScheme
          {
            viscosity_approximation,
            exact
          } peierls_creep_flow_law;

          /**
           * List of Peierls creep prefactors (A).
           */
          std::vector<double> prefactors;

          /**
           * List of Peierls creep stress exponents (n).
           */
          std::vector<double> stress_exponents;

          /**
           * List of Peierls creep activation energies (E).
           */
          std::vector<double> activation_energies;

          /**
           * List of Peierls creep activation volumes (V).
           */
          std::vector<double> activation_volumes;

          /**
           * List of Peierls stresses (sigma_p).
           */
          std::vector<double> peierls_stresses;

          /**
           * List of Peierls fitting parameters (gamma).
           */
          std::vector<double> fitting_parameters;

          /**
           * List of the first Peierls parameter related
           * to dislocation glide (p).
           */
          std::vector<double> glide_parameters_p;

          /**
           * List of the second Peierls parameter related
           * to dislocation glide (q).
           */
          std::vector<double> glide_parameters_q;

          /**
           * Parameters governing the iteration for the exact
           * Peierls viscosity.
           */
          double strain_rate_residual_threshold;
          unsigned int stress_max_iteration_number;

      };
    }
  }
}
#endif
