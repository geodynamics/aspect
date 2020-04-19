/*
  Copyright (C) 2019 by the authors of the ASPECT code.

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

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    namespace Rheology
    {
      /**
         * A data structure with all input file parameters relevant to Drucker-Prager plasticity.
         */
      struct DruckerPragerParameters
      {
        /**
         * List of static internal friction angles (phi).
         */
        std::vector<double> angles_internal_friction;

        /**
         * List of dynamic internal friction angles (phi_d).
         */
        std::vector<double> angles_internal_friction_dynamic;

        /**
         * List of the cohesions.
         */
        std::vector<double> cohesions;

        /**
         * Limit maximum yield stress from drucker prager yield criterion.
         */
        double max_yield_stress;
      };

      template <int dim>
      class DruckerPrager
      {
        public:
          /**
           * Declare the parameters this function takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters from the parameter file. These parameters might
           * be modified outside of this class for each call to the functions below,
           * so the read parameters are returned to the caller instead of stored as members.
           *
           * @param [in] n_fields The number of expected values for the angle of friction and
           * cohesion lists. Generally either the number of compositional fields or this number
           * plus one (for a background field), depending on how the user handles background fields.
           *
           * @param [in] prm The ParameterHandler to read from.
           */
          DruckerPragerParameters
          parse_parameters (const unsigned int n_fields,
                            ParameterHandler &prm);

          /**
           * Compute the steady-state angle of internal friction using the standard (static)
           * angle of internal friction, dynamic angle of internal friction, effective strain
           * rate, and reference strain rate.
           */
          double
          compute_steady_state_friction (const double angle_internal_friction,
                                         const double angle_internal_friction_dynamic,
                                         const double effective_strain_rate,
                                         const double reference_strain_rate) const;

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
                             const double max_yield_stress) const;

          /**
           * Compute the derivative of the plastic viscosity with respect to pressure.
           */
          double
          compute_derivative (const double angle_internal_friction,
                              const double effective_strain_rate) const;
      };
    }
  }
}
#endif
