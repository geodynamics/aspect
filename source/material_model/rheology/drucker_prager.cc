/*
  Copyright (C) 2019 - 2022 by the authors of the ASPECT code.

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


#include <aspect/material_model/rheology/drucker_prager.h>
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
      template <int dim>
      DruckerPrager<dim>::DruckerPrager ()
        = default;

      template <int dim>
      const DruckerPragerParameters
      DruckerPrager<dim>::compute_drucker_prager_parameters (const unsigned int composition,
                                                             const std::vector<double> &phase_function_values,
                                                             const std::vector<unsigned int> &n_phases_per_composition) const
      {
        DruckerPragerParameters drucker_prager_parameters;

        drucker_prager_parameters.max_yield_stress = max_yield_stress;

        if (phase_function_values == std::vector<double>())
          {
            // no phases
            drucker_prager_parameters.angle_internal_friction = angles_internal_friction[composition];
            drucker_prager_parameters.cohesion = cohesions[composition];
          }
        else
          {
            // Average among phases
            drucker_prager_parameters.angle_internal_friction = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phases_per_composition,
                                                                angles_internal_friction, composition);
            drucker_prager_parameters.cohesion = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phases_per_composition,
                                                 cohesions, composition);
          }
        return drucker_prager_parameters;
      }

      template <int dim>
      double
      DruckerPrager<dim>::compute_yield_stress (const double cohesion,
                                                const double angle_internal_friction,
                                                const double pressure,
                                                const double max_yield_stress) const
      {
        const double sin_phi = std::sin(angle_internal_friction);
        const double cos_phi = std::cos(angle_internal_friction);
        const double stress_inv_part = 1. / (std::sqrt(3.0) * (3.0 + sin_phi));

        // Initial yield stress (no stabilization terms)
        const double yield_stress = ( (dim==3)
                                      ?
                                      ( 6.0 * cohesion * cos_phi + 6.0 * pressure * sin_phi) * stress_inv_part
                                      :
                                      cohesion * cos_phi + pressure * sin_phi);

        return std::min(yield_stress, max_yield_stress);
      }



      template <int dim>
      double
      DruckerPrager<dim>::compute_viscosity (const double cohesion,
                                             const double angle_internal_friction,
                                             const double pressure,
                                             const double effective_strain_rate,
                                             const double max_yield_stress,
                                             const double pre_yield_viscosity) const
      {
        const double yield_stress = compute_yield_stress(cohesion, angle_internal_friction, pressure, max_yield_stress);

        const double strain_rate_effective_inv = 1./(2.*effective_strain_rate);

        double plastic_viscosity = yield_stress * strain_rate_effective_inv;

        if (use_plastic_damper)
          {
            const double total_stress = ( yield_stress + ( 2. * damper_viscosity * effective_strain_rate ) ) /
                                        ( 1 + ( damper_viscosity / pre_yield_viscosity ) );

            const double pre_yield_strain_rate = total_stress / ( 2. * pre_yield_viscosity);

            const double plastic_strain_rate = effective_strain_rate - pre_yield_strain_rate;

            plastic_viscosity = yield_stress / (2.*plastic_strain_rate) + damper_viscosity;

            // Effective viscosity
            plastic_viscosity = 1. / (1./plastic_viscosity + 1./pre_yield_viscosity);
          }

        return plastic_viscosity;
      }



      template <int dim>
      std::pair<double, double>
      DruckerPrager<dim>::compute_strain_rate_and_derivative (const double stress,
                                                              const double pressure,
                                                              const DruckerPragerParameters p) const
      {

        const double yield_stress = compute_yield_stress(p.cohesion, p.angle_internal_friction, pressure, p.max_yield_stress);

        if (stress > yield_stress)
          {
            return std::make_pair((stress - yield_stress)/(2.*damper_viscosity), 1./(2.*damper_viscosity));
          }
        else
          {
            return std::make_pair(0., 0.);
          }
      }



      template <int dim>
      double
      DruckerPrager<dim>::compute_derivative (const double angle_internal_friction,
                                              const double effective_strain_rate) const
      {
        const double sin_phi = std::sin(angle_internal_friction);

        const double stress_inv_part = 1. / (std::sqrt(3.0) * (3.0 + sin_phi));

        const double strain_rate_effective_inv = 1./(2.*effective_strain_rate);

        const double viscosity_pressure_derivative = sin_phi * strain_rate_effective_inv *
                                                     (dim == 3
                                                      ?
                                                      (6.0 * stress_inv_part)
                                                      :
                                                      1);

        return viscosity_pressure_derivative;
      }



      template <int dim>
      void
      DruckerPrager<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Angles of internal friction", "0.",
                           Patterns::Anything(),
                           "List of angles of internal friction, $\\phi$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "For a value of zero, in 2D the von Mises criterion is retrieved. "
                           "Angles higher than 30 degrees are harder to solve numerically. Units: degrees.");
        prm.declare_entry ("Cohesions", "1e20",
                           Patterns::Anything(),
                           "List of cohesions, $C$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "The extremely large default cohesion value (1e20 Pa) prevents the viscous stress from "
                           "exceeding the yield stress. Units: \\si{\\pascal}.");
        prm.declare_entry ("Maximum yield stress", "1e12", Patterns::Double (0.),
                           "Limits the maximum value of the yield stress determined by the "
                           "Drucker-Prager plasticity parameters. Default value is chosen so this "
                           "is not automatically used. Values of 100e6--1000e6 $Pa$ have been used "
                           "in previous models. Units: \\si{\\pascal}.");
        prm.declare_entry ("Use plastic damper","false",
                           Patterns::Bool (),
                           "Whether to use a plastic damper when computing the Drucker-Prager "
                           "plastic viscosity. The damper acts to stabilize the plastic shear "
                           "band width and remove associated mesh-dependent behavior at "
                           "sufficient resolutions.");
        prm.declare_entry ("Plastic damper viscosity", "0.0", Patterns::Double(0),
                           "Viscosity of the damper that acts in parallel with the plastic viscosity "
                           "to produce mesh-independent behavior at sufficient resolutions. Units: \\si{\\pascal\\second}");
      }



      template <int dim>
      void
      DruckerPrager<dim>::parse_parameters (ParameterHandler &prm,
                                            const std::unique_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition)
      {
        // Retrieve the list of composition names
        const std::vector<std::string> list_of_composition_names = this->introspection().get_composition_names();
        // Establish that a background field is required here
        const bool has_background_field = true;

        angles_internal_friction = Utilities::parse_map_to_double_array(prm.get("Angles of internal friction"),
                                                                        list_of_composition_names,
                                                                        has_background_field,
                                                                        "Angles of internal friction",
                                                                        true,
                                                                        expected_n_phases_per_composition);

        // Convert angles from degrees to radians
        for (double &angle : angles_internal_friction)
          angle *= numbers::PI/180.0;

        cohesions = Utilities::parse_map_to_double_array(prm.get("Cohesions"),
                                                         list_of_composition_names,
                                                         has_background_field,
                                                         "Cohesions",
                                                         true,
                                                         expected_n_phases_per_composition);

        // Limit maximum value of the Drucker-Prager yield stress
        max_yield_stress = prm.get_double("Maximum yield stress");

        // Whether to include a plastic damper when computing the Drucker-Prager plastic viscosity
        use_plastic_damper = prm.get_bool("Use plastic damper");

        // Stabilize plasticity through a viscous damper.
        // The viscosity of the damper is implicitly zero if it is not used
        if (use_plastic_damper)
          damper_viscosity = prm.get_double("Plastic damper viscosity");
        else
          damper_viscosity = 0.;

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
    template class DruckerPrager<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
