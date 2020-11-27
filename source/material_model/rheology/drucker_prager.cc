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


#include <aspect/material_model/rheology/drucker_prager.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
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
        double yield_stress = ( (dim==3)
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

        if (use_plastic_damper == true)
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
                           Patterns::List(Patterns::Double (0.)),
                           "List of angles of internal friction, $\\phi$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "For a value of zero, in 2D the von Mises criterion is retrieved. "
                           "Angles higher than 30 degrees are harder to solve numerically. Units: degrees.");
        prm.declare_entry ("Cohesions", "1e20",
                           Patterns::List(Patterns::Double (0.)),
                           "List of cohesions, $C$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "The extremely large default cohesion value (1e20 Pa) prevents the viscous stress from "
                           "exceeding the yield stress. Units: \\si{\\pascal}.");
        prm.declare_entry ("Maximum yield stress", "1e12", Patterns::Double (0.),
                           "Limits the maximum value of the yield stress determined by the "
                           "drucker-prager plasticity parameters. Default value is chosen so this "
                           "is not automatically used. Values of 100e6--1000e6 $Pa$ have been used "
                           "in previous models. Units: \\si{\\pascal}.");
        prm.declare_entry ("Use plastic damper","false",
                           Patterns::Bool (),
                           "Whether to use a plastic damper when computing the drucker-prager "
                           "plastic viscosity. The damper acts to stabilize the plastic shear "
                           "band width and remove associated mesh-dependent behavior at "
                           "sufficient resolutions.");
        prm.declare_entry ("Plastic damper viscosity", "0.0", Patterns::Double(0),
                           "Viscous damper that acts in parallel with the plastic viscosity "
                           "to produce mesh-independent behavior at sufficient resolutions. Units: \\si{\\pascal\\second}");
      }



      template <int dim>
      DruckerPragerParameters
      DruckerPrager<dim>::parse_parameters (const unsigned int n_fields,
                                            ParameterHandler &prm)
      {
        DruckerPragerParameters parameters;

        parameters.angles_internal_friction = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Angles of internal friction"))),
                                                                                      n_fields,
                                                                                      "Angles of internal friction");
        // Convert angles from degrees to radians
        for (double &angle : parameters.angles_internal_friction)
          angle *= numbers::PI/180.0;

        parameters.cohesions = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Cohesions"))),
                                                                       n_fields,
                                                                       "Cohesions");

        // Limit maximum value of the drucker-prager yield stress
        parameters.max_yield_stress = prm.get_double("Maximum yield stress");

        // Whether to include a plastic damper when computing the drucker-prager plastic viscosity
        use_plastic_damper = prm.get_bool("Use plastic damper");

        // Stabalize plasticity through a viscous damper
        damper_viscosity = prm.get_double("Plastic damper viscosity");

        return parameters;
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
