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


#include <aspect/material_model/rheology/diffusion_creep.h>
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
      DiffusionCreep<dim>::DiffusionCreep ()
      {}



      template <int dim>
      double
      DiffusionCreep<dim>::compute_viscosity (const double pressure,
                                              const double temperature,
                                              const unsigned int composition) const
      {
        // Power law creep equation
        //    viscosity = 0.5 * A^(-1/n) * d^(m/n) * exp((E + P*V)/(nRT))
        // A: prefactor,
        // d: grain size, m: grain size exponent, E: activation energy, P: pressure,
        // V; activation volume, R: gas constant, T: temperature.
        const double viscosity_diffusion = 0.5 / prefactors_diffusion[composition] *
                                           std::exp((activation_energies_diffusion[composition] +
                                                     pressure*activation_volumes_diffusion[composition])/
                                                    (constants::gas_constant*temperature)) *
                                           std::pow(grain_size, grain_size_exponents_diffusion[composition]);

        return viscosity_diffusion;
      }



      template <int dim>
      void
      DiffusionCreep<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Prefactors for diffusion creep", "1.5e-15",
                           Patterns::List(Patterns::Double (0.)),
                           "List of viscosity prefactors, $A$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value. "
                           "Units: $Pa^{-1} m^{m_{\\text{diffusion}}} s^{-1}$");
        prm.declare_entry ("Grain size exponents for diffusion creep", "3.",
                           Patterns::List(Patterns::Double (0.)),
                           "List of grain size exponents, $m_{\\text{diffusion}}$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value. Units: None");
        prm.declare_entry ("Activation energies for diffusion creep", "375e3",
                           Patterns::List(Patterns::Double (0.)),
                           "List of activation energies, $E_a$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value.  Units: $J / mol$");
        prm.declare_entry ("Activation volumes for diffusion creep", "6e-6",
                           Patterns::List(Patterns::Double (0.)),
                           "List of activation volumes, $V_a$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value.  Units: $m^3 / mol$");
        prm.declare_entry ("Grain size", "1e-3", Patterns::Double (0.), "Units: $m$");
      }



      template <int dim>
      void
      DiffusionCreep<dim>::parse_parameters (ParameterHandler &prm)
      {
        // increment by one for background:
        const unsigned int n_fields = this->n_compositional_fields() + 1;

        prefactors_diffusion = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Prefactors for diffusion creep"))),
                                                                       n_fields,
                                                                       "Prefactors for diffusion creep");
        grain_size_exponents_diffusion = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Grain size exponents for diffusion creep"))),
                                                                                 n_fields,
                                                                                 "Grain size exponents for diffusion creep");
        activation_energies_diffusion = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Activation energies for diffusion creep"))),
                                                                                n_fields,
                                                                                "Activation energies for diffusion creep");
        activation_volumes_diffusion = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Activation volumes for diffusion creep"))),
                                                                               n_fields,
                                                                               "Activation volumes for diffusion creep");
        grain_size = prm.get_double("Grain size");
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
    template class DiffusionCreep<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}
