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


#include <aspect/material_model/rheology/drucker_prager_power.h>
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
      DruckerPragerPower<dim>::DruckerPragerPower ()
        = default;



      template <int dim>
      const DruckerPragerParameters
      DruckerPragerPower<dim>::compute_drucker_prager_parameters (const unsigned int composition,
                                                                  const std::vector<double> &phase_function_values,
                                                                  const std::vector<unsigned int> &n_phase_transitions_per_composition) const
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
            drucker_prager_parameters.angle_internal_friction = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                                                angles_internal_friction, composition);
            drucker_prager_parameters.cohesion = MaterialModel::MaterialUtilities::phase_average_value(phase_function_values, n_phase_transitions_per_composition,
                                                 cohesions, composition);
          }
        return drucker_prager_parameters;
      }



      template <int dim>
      double
      DruckerPragerPower<dim>::compute_yield_stress (const double cohesion,
                                                     const double angle_internal_friction,
                                                     const double pressure,
                                                     const double max_yield_stress) const
      {
        const double sin_phi = std::sin(angle_internal_friction);
        const double cos_phi = std::cos(angle_internal_friction);

        // The expression below differs from Eq. 9 of Glerum et al, 2018.
        // There are actually three different ways of choosing this parameter, which
        // correspond to the Drucker-Prager yield surface either
        // circumscribing (Glerum et al 2018), middle circumscribing or
        // inscribing the Mohr-Coulomb yield surface.
        // See for instance Owen & Hinton, Finite Elements in Plasticity, 1980.
        // Here the middle circumscribing approach is taken.
        const double stress_inv_part = 1. / (std::sqrt(3.0) * (3.0 + sin_phi));

        // Initial yield stress (no stabilization terms)
        const double yield_stress = ( (dim==3)
                                      ?
                                      (( 6.0 * cohesion * cos_phi + 6.0 * pressure * sin_phi) * stress_inv_part)
                                      :
                                      (cohesion * cos_phi + pressure * sin_phi));

        return std::min(yield_stress, max_yield_stress);
      }



      template <int dim>
      double
      DruckerPragerPower<dim>::compute_viscosity (const double cohesion,
                                                  const double angle_internal_friction,
                                                  const double pressure,
                                                  const double effective_strain_rate,
                                                  const double max_yield_stress) const
      {
        const double yield_stress = compute_yield_stress(cohesion, angle_internal_friction, pressure, max_yield_stress);
        const double stress = yield_stress * std::pow(effective_strain_rate/drucker_prager_edot_ref, 1./drucker_prager_stress_exponent);
        const double apparent_viscosity = stress / (2. * effective_strain_rate);

        return apparent_viscosity;
      }



      template <int dim>
      std::pair<double, double>
      DruckerPragerPower<dim>::compute_strain_rate_and_derivative (const double stress,
                                                                   const double pressure,
                                                                   const DruckerPragerParameters &p) const
      {

        const double yield_stress = compute_yield_stress(p.cohesion, p.angle_internal_friction, pressure, p.max_yield_stress);
        const double strain_rate = drucker_prager_edot_ref * std::pow(stress/yield_stress, drucker_prager_stress_exponent);
        const double deriv = drucker_prager_stress_exponent * strain_rate / stress;
        return std::make_pair(strain_rate, deriv);
      }



      template <int dim>
      std::pair<double, double>
      DruckerPragerPower<dim>::compute_log_strain_rate_and_derivative (const double log_stress,
                                                                       const double pressure,
                                                                       const DruckerPragerParameters &p) const
      {

        const double yield_stress = compute_yield_stress(p.cohesion, p.angle_internal_friction, pressure, p.max_yield_stress);
        const double log_strain_rate = drucker_prager_log_edot_ref + drucker_prager_stress_exponent * (log_stress - std::log(yield_stress));

        return std::make_pair(log_strain_rate, drucker_prager_stress_exponent);

      }



      template <int dim>
      void
      DruckerPragerPower<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Angles of internal friction", "0.",
                           Patterns::Anything(),
                           "List of angles of internal friction, $\\phi$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "For a value of zero, in 2d the von Mises criterion is retrieved. "
                           "Angles higher than 30 degrees are harder to solve numerically. Units: degrees.");
        prm.declare_entry ("Cohesions", "1e20",
                           Patterns::Anything(),
                           "List of cohesions, $C$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "The extremely large default cohesion value (1e20 Pa) prevents the viscous stress from "
                           "exceeding the yield stress. Units: \\si{\\pascal}.");
        prm.declare_entry ("Maximum yield stress", "1e12", Patterns::Double (0.),
                           "Limits the maximum value of the yield stress determined by the "
                           "Drucker-Prager plasticity parameters. Default value is chosen so this "
                           "is not automatically used. Values of 100e6--1000e6 $Pa$ have been used "
                           "in previous models. Units: \\si{\\pascal}.");
        prm.declare_entry ("Reference plastic strain rate", "1e-18", Patterns::Double (0.),
                           "Provides the strain rate at which the yield stress determined by the "
                           "Drucker-Prager plasticity parameters apply. Units: \\si{\\per\\second}.");
        prm.declare_entry ("Plastic stress exponent", "50", Patterns::Double (0.),
                           "Provides the stress exponent that modifies the yield stress according to "
                           "the strain rate. The default value is chosen to provide trade-off between "
                           "a sharp onset of plasticity and smooth transition between flow mechanisms "
                           "required for iterative strain rate decomposition. Units: None.");
      }



      template <int dim>
      void
      DruckerPragerPower<dim>::parse_parameters (ParameterHandler &prm,
                                                 const std::unique_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition)
      {
        // Retrieve the list of composition names
        std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();

        // Retrieve the list of names of fields that represent chemical compositions, and not, e.g.,
        // plastic strain
        std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();

        // Establish that a background field is required here
        compositional_field_names.insert(compositional_field_names.begin(), "background");
        chemical_field_names.insert(chemical_field_names.begin(), "background");

        Utilities::MapParsing::Options options(chemical_field_names, "Angles of internal friction");
        options.list_of_allowed_keys = compositional_field_names;

        if (expected_n_phases_per_composition)
          {
            options.allow_multiple_values_per_key = true;
            options.n_values_per_key = *expected_n_phases_per_composition;

            // check_values_per_key is required to be true to duplicate single values
            // if they are to be used for all phases associated with a given key.
            options.check_values_per_key = true;
          }

        angles_internal_friction = Utilities::MapParsing::parse_map_to_double_array(prm.get("Angles of internal friction"),
                                                                                    options);

        // Convert angles from degrees to radians
        for (double &angle : angles_internal_friction)
          angle *= constants::degree_to_radians;

        options.property_name = "Cohesions";
        cohesions = Utilities::MapParsing::parse_map_to_double_array(prm.get("Cohesions"),
                                                                     options);

        // Limit maximum value of the Drucker-Prager yield stress
        max_yield_stress = prm.get_double("Maximum yield stress");

        drucker_prager_edot_ref = prm.get_double("Reference plastic strain rate");
        drucker_prager_log_edot_ref = std::log(drucker_prager_edot_ref);
        drucker_prager_stress_exponent = prm.get_double("Plastic stress exponent");

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
    template class DruckerPragerPower<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
