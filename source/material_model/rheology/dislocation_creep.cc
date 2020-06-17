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


#include <aspect/material_model/rheology/dislocation_creep.h>
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
      DislocationCreep<dim>::DislocationCreep ()
      {}



      template <int dim>
      double
      DislocationCreep<dim>::compute_viscosity (const double strain_rate,
                                                const double pressure,
                                                const double temperature,
                                                const unsigned int composition) const
      {
        // Power law creep equation:
        //    viscosity = 0.5 * A^(-1/n) * edot_ii^((1-n)/n) * exp((E + P*V)/(nRT))
        // A: prefactor, edot_ii: square root of second invariant of deviatoric strain rate tensor,
        // E: activation energy, P: pressure,
        // V; activation volume, n: stress exponent, R: gas constant, T: temperature.
        const double viscosity_dislocation = 0.5 * std::pow(prefactors_dislocation[composition],-1/stress_exponents_dislocation[composition]) *
                                             std::exp((activation_energies_dislocation[composition] + pressure*activation_volumes_dislocation[composition])/
                                                      (constants::gas_constant*temperature*stress_exponents_dislocation[composition])) *
                                             std::pow(strain_rate,((1. - stress_exponents_dislocation[composition])/stress_exponents_dislocation[composition]));
        return viscosity_dislocation;
      }



      template <int dim>
      void
      DislocationCreep<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Prefactors for dislocation creep", "1.1e-16",
                           Patterns::List(Patterns::Double (0.)),
                           "List of viscosity prefactors, $A$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value. "
                           "Units: $Pa^{-n_{\\text{dislocation}}} s^{-1}$");
        prm.declare_entry ("Stress exponents for dislocation creep", "3.5",
                           Patterns::List(Patterns::Double (0.)),
                           "List of stress exponents, $n_{\\text{dislocation}}$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value.  Units: None");
        prm.declare_entry ("Activation energies for dislocation creep", "530e3",
                           Patterns::List(Patterns::Double (0.)),
                           "List of activation energies, $E_a$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value.  Units: $J / mol$");
        prm.declare_entry ("Activation volumes for dislocation creep", "1.4e-5",
                           Patterns::List(Patterns::Double (0.)),
                           "List of activation volumes, $V_a$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value.  Units: $m^3 / mol$");
      }



      template <int dim>
      void
      DislocationCreep<dim>::parse_parameters (ParameterHandler &prm)
      {
        // increment by one for background:
        const unsigned int n_fields = this->n_compositional_fields() + 1;

        prefactors_dislocation = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Prefactors for dislocation creep"))),
                                                                         n_fields,
                                                                         "Prefactors for dislocation creep");
        stress_exponents_dislocation = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Stress exponents for dislocation creep"))),
                                                                               n_fields,
                                                                               "Stress exponents for dislocation creep");
        activation_energies_dislocation = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Activation energies for dislocation creep"))),
                                                                                  n_fields,
                                                                                  "Activation energies for dislocation creep");
        activation_volumes_dislocation = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Activation volumes for dislocation creep"))),
                                                                                 n_fields,
                                                                                 "Activation volumes for dislocation creep");
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
    template class DislocationCreep<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}
