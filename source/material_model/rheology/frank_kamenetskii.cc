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


#include <aspect/material_model/rheology/frank_kamenetskii.h>
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
      FrankKamenetskii<dim>::FrankKamenetskii ()
      {}



      template <int dim>
      double
      FrankKamenetskii<dim>::compute_viscosity (const double temperature,
                                                const unsigned int composition) const
      {
        const double reference_temperature = this->get_adiabatic_surface_temperature();

        const double viscosity_frank_kamenetskii = prefactors_frank_kamenetskii[composition] * std::exp(viscosity_ratios_frank_kamenetskii[composition] * 0.5 * (1.0-temperature/reference_temperature));

        return viscosity_frank_kamenetskii;
      }





      template <int dim>
      void
      FrankKamenetskii<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Viscosity ratios for Frank Kamenetskii", "15.",
                           Patterns::List(Patterns::Double (0.)),
                           "An adjusted viscosity ratio, $E$, for the viscosity approximation, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value. "
                           "Units: None");
        prm.declare_entry ("Prefactors for Frank Kamenetskii", "1.e21",
                           Patterns::List(Patterns::Double (0.)),
                           "A viscosity prefactor for the viscosity approximation, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value.  Units: None");
      }



      template <int dim>
      void
      FrankKamenetskii<dim>::parse_parameters (ParameterHandler &prm)
      {
        // increment by one for background:
        const unsigned int n_fields = this->n_compositional_fields() + 1;

        AssertThrow (this->include_adiabatic_heating() == false,
                     ExcMessage("The Frank-Kamenetskii rheology is currently only implemented for "
                                "models without adiabatic heating. Please implement the necessary "
                                "temperature adjustment if you need this feature."));

        AssertThrow (this->get_adiabatic_surface_temperature() > 0.0,
                     ExcMessage("The Frank-Kamenetskii rheology can only be used when the adiabatic "
                                "surface temperature (reference_temperature in equation for viscosity) "
                                "is non-zero."));

        viscosity_ratios_frank_kamenetskii = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Viscosity ratios for Frank Kamenetskii"))),
                                                                                     n_fields,
                                                                                     "Viscosity ratios for Frank Kamenetskii");
        prefactors_frank_kamenetskii = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Prefactors for Frank Kamenetskii"))),
                                                                               n_fields,
                                                                               "Prefactors for Frank Kamenetskii");
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
    template class FrankKamenetskii<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
