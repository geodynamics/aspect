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


#include <aspect/material_model/rheology/drucker_prager.h>
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
      {}



      template <int dim>
      double
      DruckerPrager<dim>::compute_yield_stress (const double cohesion,
                                                const double angle_internal_friction,
                                                const double pressure) const
      {
        const double sin_phi = std::sin(angle_internal_friction);
        const double cos_phi = std::cos(angle_internal_friction);
        const double stress_inv_part = 1. / (std::sqrt(3.0) * (3.0 + sin_phi));

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
                                             const double effective_strain_rate) const
      {
        const double yield_stress = compute_yield_stress(cohesion, angle_internal_friction, pressure);

        const double strain_rate_effective_inv = 1./(2.*effective_strain_rate);

        return yield_stress * strain_rate_effective_inv;
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
      double
      DruckerPrager<dim>::get_max_yield_stress () const
      {
        return max_yield_stress;
      }



      template <int dim>
      std::vector<double>
      DruckerPrager<dim>::get_cohesions () const
      {
        return cohesions;
      }



      template <int dim>
      std::vector<double>
      DruckerPrager<dim>::get_angles_internal_friction () const
      {
        return angles_internal_friction;
      }


      template <int dim>
      void
      DruckerPrager<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Angles of internal friction", "0",
                           Patterns::List(Patterns::Double(0)),
                           "List of angles of internal friction, $\\phi$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "For a value of zero, in 2D the von Mises criterion is retrieved. "
                           "Angles higher than 30 degrees are harder to solve numerically. Units: degrees.");
        prm.declare_entry ("Cohesions", "1e20",
                           Patterns::List(Patterns::Double(0)),
                           "List of cohesions, $C$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "The extremely large default cohesion value (1e20 Pa) prevents the viscous stress from "
                           "exceeding the yield stress. Units: $Pa$.");
        prm.declare_entry ("Maximum yield stress", "1e12", Patterns::Double(0),
                           "Limits the maximum value of the yield stress determined by the "
                           "drucker-prager plasticity parameters. Default value is chosen so this "
                           "is not automatically used. Values of 100e6--1000e6 $Pa$ have been used "
                           "in previous models. Units: $Pa$");
      }



      template <int dim>
      void
      DruckerPrager<dim>::parse_parameters (ParameterHandler &prm)
      {
        // increment by one for background:
        const unsigned int n_fields = this->n_compositional_fields() + 1;

        angles_internal_friction = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Angles of internal friction"))),
                                                                           n_fields,
                                                                           "Angles of internal friction");
        // Convert angles from degrees to radians
        for (unsigned int i = 0; i<n_fields; ++i)
          angles_internal_friction[i] *= numbers::PI/180.0;

        cohesions = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Cohesions"))),
                                                            n_fields,
                                                            "Cohesions");

        // Limit maximum value of the drucker-prager yield stress
        max_yield_stress = prm.get_double("Maximum yield stress");
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
  }
}
