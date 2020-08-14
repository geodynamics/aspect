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


#include <aspect/material_model/rheology/peierls_creep.h>
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
      PeierlsCreep<dim>::PeierlsCreep ()
      {}



      template <int dim>
      double
      PeierlsCreep<dim>::compute_viscosity (const double strain_rate,
                                            const double pressure,
                                            const double temperature,
                                            const unsigned int composition) const
      {
        const unsigned int c = composition;

        double viscosity = 0.0;

        switch (peierls_creep_flow_law)
          {
            case viscosity_approximation:
            {
              /**
               * An approximation of the Peierls creep formulation, where stress is replaced with strain rate
               * (second invariant) when calculating viscosity. This substitution requires a fitting parameter (gamma),
               * which describes the relationship between stress and the Peierls stress (gamma = stress/stress_peierls).
               * Details of this derivation can be found in the following document:
               * https://ucdavis.app.box.com/s/cl5mwhkjeabol4otrdukfcdwfvg9me4w/file/705438695737
               * The formulation for the viscosity is:
               *   stress_term * arrhenius_term * strain_rate_term, where
               * stress_term      = (0.5 * gamma * sigma_p) / (A * ((gamma * sigma_p)^n)^(1/(s+n)))
               * arrhenius_term   = exp(( (E + P * V) / (R * T)) * (1 - gamma^p)^q / (s + n) )
               * strain_rate_term = edot_ii^(1 / (s + n) - 1)
               * Above,
               * s = ((E + P * V) / (R * T)) * p * q * (1 - gamma^p)^(q-1) * (gamma^p)
               * sigma_p is the Peierls stress
               * gamma is the Peierls creep fitting parameter
               * p is the first Peierls glide parameter
               * q is the second Peierls creep glide parameter
               * A is the Peierls prefactor term (1/s 1/Pa^2)
               * n is the Peierls stress exponent
               * E is the Peierls activation energy (J/mol)
               * V is the Peierls activation volume (m^3/mol). However, to date V has always been set to 0.
               * edot_ii is the second invariant of the deviatoric strain rate (1/s)
               * T is temperature (K)
               * R is the gas constant
               */
              const double s = ( (activation_energies[c] + pressure * activation_volumes[c]) / (constants::gas_constant * temperature)) *
                               glide_parameters_p[c] * glide_parameters_q[c] *
                               std::pow((1. - std::pow(fitting_parameters[c], glide_parameters_p[c])),(glide_parameters_q[c] - 1.)) *
                               std::pow(fitting_parameters[c], glide_parameters_p[c]);

              const double stress_term = 0.5 * (fitting_parameters[c] * peierls_stresses[c]) /
                                         std::pow((prefactors[c] * std::pow(fitting_parameters[c] * peierls_stresses[c],stress_exponents[c])),( 1. / (s + stress_exponents[c])));

              const double arrhenius_term = std::exp( ((activation_energies[c] + pressure * activation_volumes[c]) / (constants::gas_constant * temperature)) *
                                                      (std::pow((1. - std::pow(fitting_parameters[c],glide_parameters_p[c])),glide_parameters_q[c])) /
                                                      (s + stress_exponents[c]) );

              const double strain_rate_term = std::pow(strain_rate, (1. / (s + stress_exponents[c])) - 1.);

              viscosity = stress_term * arrhenius_term * strain_rate_term;
              break;
            }
            default:
            {
              AssertThrow(false, ExcNotImplemented());
              break;
            }
          }

        return viscosity;
      }



      template <int dim>
      void
      PeierlsCreep<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Peierls creep flow law", "viscosity approximation",
                           Patterns::Selection("viscosity approximation"),
                           "Select what type of Peierls creep flow law to use. Currently, the "
                           "only available option is an approximation to Peierls creep, which "
                           "uses the strain rate invariant, rather than stress, as an input. ");
        prm.declare_entry ("Prefactors for Peierls creep", "1.4e-19",
                           Patterns::List(Patterns::Double(0.)),
                           "List of viscosity prefactors, $A$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\pascal}$^{-n_{\\text{peierls}}}$ \\si{\\per\\second}");
        prm.declare_entry ("Stress exponents for Peierls creep", "2.0",
                           Patterns::List(Patterns::Double(0.)),
                           "List of stress exponents, $n_{\\text{peierls}}$, for background material and compositional "
                           "fields, for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value.  Units: None.");
        prm.declare_entry ("Activation energies for Peierls creep", "320e3",
                           Patterns::List(Patterns::Double(0)),
                           "List of activation energies, $E$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value. Units: \\si{\\joule\\per\\mole}.");
        prm.declare_entry ("Activation volumes for Peierls creep", "1.4e-5",
                           Patterns::List(Patterns::Double(0.)),
                           "List of activation volumes, $V$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\meter\\cubed\\per\\mole}.");
        prm.declare_entry ("Peierls stresses", "5.e9",
                           Patterns::List(Patterns::Double(0.)),
                           "List of stress limits for Peierls creep $\\sigma_{\\text{peierls}}$ for background "
                           "material and compositional fields, for a total of N+1 values, where N is the number "
                           "of compositional fields. If only one value is given, then all use the same value. "
                           "Units: \\si{\\pascal}");
        prm.declare_entry ("Peierls fitting parameters", "0.17",
                           Patterns::List(Patterns::Double(0.)),
                           "List of fitting parameters $\\gamma$ between stress $\\sigma$ and the Peierls "
                           "stress $\\sigma_{\\text{peierls}}$ for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. If only one "
                           "value is given, then all use the same value. Units: none");
        prm.declare_entry ("Peierls glide parameters p", "0.5",
                           Patterns::List(Patterns::Double(0.)),
                           "List of the first Peierls creep glide parameters, $p$, for background and compositional "
                           "fields for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value. Units: none");
        prm.declare_entry ("Peierls glide parameters q", "1.0",
                           Patterns::List(Patterns::Double(0.)),
                           "List of the second Peierls creep glide parameters, $q$, for background and compositional "
                           "fields for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value. Units: none");

      }



      template <int dim>
      void
      PeierlsCreep<dim>::parse_parameters (ParameterHandler &prm)
      {
        // increment by one for background:
        const unsigned int n_fields = this->n_compositional_fields() + 1;

        if (prm.get ("Peierls creep flow law") == "viscosity approximation")
          peierls_creep_flow_law = viscosity_approximation;
        else
          AssertThrow(false, ExcMessage("Not a valid Peierls creep flow law"));

        prefactors = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Prefactors for Peierls creep"))),
                                                             n_fields,
                                                             "Prefactors for Peierls creep");

        stress_exponents = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Stress exponents for Peierls creep"))),
                                                                   n_fields,
                                                                   "Stress exponents for Peierls creep");

        activation_energies = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Activation energies for Peierls creep"))),
                                                                      n_fields,
                                                                      "Activation energies for Peierls creep");

        activation_volumes = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Activation volumes for Peierls creep"))),
                                                                     n_fields,
                                                                     "Activation volumes for Peierls creep");

        peierls_stresses = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Peierls stresses"))),
                                                                   n_fields,
                                                                   "Peierls stresses");

        fitting_parameters = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Peierls fitting parameters"))),
                                                                     n_fields,
                                                                     "Peierls fitting parameters");

        glide_parameters_p = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Peierls glide parameters p"))),
                                                                     n_fields,
                                                                     "Peierls glide parameters p");

        glide_parameters_q = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Peierls glide parameters q"))),
                                                                     n_fields,
                                                                     "Peierls glide parameters q");
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
    template class PeierlsCreep<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
