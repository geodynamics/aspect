/*
  Copyright (C) 2019 - 2020 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_rheology_drucker_prager_h
#define _aspect_material_model_rheology_drucker_prager_h

#include <aspect/global.h>
#include <aspect/material_model/interface.h>

#include <aspect/material_model/utilities.h>
#include <aspect/simulator_access.h>
namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    namespace Rheology
    {
      /**
       * A data structure with all input file parameters relevant to
       * Drucker-Prager plasticity.
       */
      struct DruckerPragerParameters
      {
        /**
         * Internal friction angle (phi) for the current composition and phase
         */
        double angle_internal_friction;

        /**
         * Cohesion for the current composition and phase
         */
        double cohesion;

        /**
         * Limit maximum yield stress from drucker prager yield criterion.
         */
        double max_yield_stress;
      };

      template <int dim>
      class DruckerPrager : public ::aspect::SimulatorAccess<dim>
      {
        public:
          DruckerPrager();
          /**
           * Declare the parameters this function takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           * If @p expected_n_phases_per_composition points to a vector of
           * unsigned integers, this is considered the number of phases
           * for each compositional field (plus possibly a background field)
           * and this number will be checked against the parsed parameters.
           *
           * @param [in] prm The ParameterHandler to read from.
           */
          void
          parse_parameters (ParameterHandler &prm,
                            const std::unique_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition = nullptr);

          const DruckerPragerParameters
          compute_drucker_prager_parameters (const unsigned int composition,
                                             const std::vector<double> &phase_function_values = std::vector<double>(),
                                             const std::vector<unsigned int> &n_phases_per_composition = std::vector<unsigned int>()) const;

          /**
           * Compute the plastic yield stress based on the Drucker Prager yield criterion.
           */
          double
          compute_yield_stress (const double cohesion,
                                const double angle_internal_friction,
                                const double pressure,
                                const double max_yield_stress) const;

          /**
           * Compute the plastic viscosity with the yield stress and effective strain rate.
           */
          double
          compute_viscosity (const double cohesion,
                             const double angle_internal_friction,
                             const double pressure,
                             const double effective_strain_rate,
                             const double max_yield_stress,
                             const double pre_yield_viscosity = std::numeric_limits<double>::infinity()) const;

          /**
           * Compute the strain rate and first stress derivative
           * as a function of stress based on the damped Drucker-Prager flow law.
           */
          std::pair<double, double>
          compute_strain_rate_and_derivative (const double stress,
                                              const double pressure,
                                              const DruckerPragerParameters p) const;

          /**
           * Compute the derivative of the plastic viscosity with respect to pressure.
           */
          double
          compute_derivative (const double angle_internal_friction,
                              const double effective_strain_rate) const;

        private:

          std::vector<double> angles_internal_friction;
          std::vector<double> cohesions;
          double max_yield_stress;

          /**
           * Whether to add a plastic damper in the computation
           * of the drucker prager plastic viscosity.
           */
          bool use_plastic_damper;

          /**
           * Viscosity of a damper used to stabilize plasticity
           */
          double damper_viscosity;
      };
    }
  }
}
#endif
