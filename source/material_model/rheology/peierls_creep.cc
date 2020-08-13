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
            case viscosity_approx:
            {
              /**
               * An approximation of the Peierls creep formulation, where stress is replaced with strain rate
               * (second invariant) when calculating viscosity. This substitution requires three additional
               * fitting parameters (gamma,p,q). Details of this derivation can be found in the following document:
               * https://ucdavis.app.box.com/s/cl5mwhkjeabol4otrdukfcdwfvg9me4w/file/705438695737
               * The formulation for the viscosity is:
               *  ((0.5*gamma*sigma_p)/(A*((gamma*sigma_p)^n)^(1/(s+n)))) * exp(((E+PV)/(R*T))*(((1-gamma^p)^q)/(s+n))) * edot_ii^(1/(s+n)-1)
               * where
               * s = ((E + P*V)/(R*T))*p*q*((1 - gamma^p)^(q-1)*(gamma^p)
               * sigma_p is the Peierls stress
               * gamma is the Peierls creep fitting parameter
               * p is the first Peierls creep fitting exponent
               * q is the second Peierls creep fitting exponent
               * A is the Peierls prefactor term (1/s 1/Pa^2)
               * n is the Peierls stress exponent
               * E is the Peierls activation energy (J/mol)
               * V is the Peierls activation volume (m^3/mol)
               * edot_ii is the second invariant of the deviatoric strain rate (1/s)
               * T is temperature (K)
               * R is the gas constant
               */
              const double s = ( (activation_energies_peierls[c] + pressure * activation_volumes_peierls[c]) / (constants::gas_constant * temperature)) *
                               peierls_fitting_exponents_p[c] * peierls_fitting_exponents_q[c] *
                               std::pow((1. - std::pow(peierls_fitting_parameters[c], peierls_fitting_exponents_p[c])),(peierls_fitting_exponents_q[c] - 1.0)) *
                               std::pow(peierls_fitting_parameters[c], peierls_fitting_exponents_p[c]);

              const double peierls_left_term = 0.5 * (peierls_fitting_parameters[c] * peierls_stresses[c]) /
                                               std::pow((prefactors_peierls[c] * std::pow(peierls_fitting_parameters[c] * peierls_stresses[c],stress_exponents_peierls[c])),( 1 / (s + stress_exponents_peierls[c])));

              const double peierls_middle_term_a = (activation_energies_peierls[c] + pressure * activation_volumes_peierls[c]) / (constants::gas_constant * temperature);

              const double peierls_middle_term_b = std::pow((1. - std::pow(peierls_fitting_parameters[c],peierls_fitting_exponents_p[c])),peierls_fitting_exponents_q[c]);

              const double peierls_middle_term_c = s + stress_exponents_peierls[c];

              const double peierls_middle_term = std::exp(peierls_middle_term_a * peierls_middle_term_b / peierls_middle_term_c);

              const double peierls_right_term = std::pow(strain_rate, (1. / (s + stress_exponents_peierls[c])) - 1.);

              viscosity = peierls_left_term * peierls_middle_term * peierls_right_term;
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
                           "uses the strain rate invariant, rather than stress, as an input. "
                           "Future options will allow formulations that use stress as an input.");
        prm.declare_entry ("Prefactors for Peierls creep", "1.4e-7",
                           Patterns::List(Patterns::Double(0)),
                           "List of viscosity prefactors, $A$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value. Units: $1/s 1/Pa^{2}$");
        prm.declare_entry ("Stress exponents for Peierls creep", "2.0",
                           Patterns::Anything(),
                           "List of stress exponents, $n_{\\text{peierls}}$, for background material and compositional "
                           "fields, for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value.  Units: None.");
        prm.declare_entry ("Activation energies for Peierls creep", "320e3",
                           Patterns::List(Patterns::Double(0)),
                           "List of activation energies, $E$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value.  Units: $J / mol$");
        prm.declare_entry ("Activation volumes for Peierls creep", "1.4e-5",
                           Patterns::Anything(),
                           "List of activation volumes, $V$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value. "
                           "Units: \\si{\\meter\\cubed\\per\\mole}.");
        prm.declare_entry ("Peierls stresses", "5.e9",
                           Patterns::List(Patterns::Double(0)),
                           "List of stress limits for Peierls creep $\\sigma_{\\text{peierls}}$ for background "
                           "material and compositional fields, for a total of N+1 values, where N is the number "
                           "of compositional fields. If only one value is given, then all use the same value. "
                           "Units: $Pa$");
        prm.declare_entry ("Peierls fitting parameters", "0.17",
                           Patterns::List(Patterns::Double(0)),
                           "List of fitting parameters $\\gamma$ between stress $\\sigma$ and the Peierls "
                           "stress $\\sigma_{\\text{peierls}}$ for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. If only one "
                           "value is given, then all use the same value. Units: none");
        prm.declare_entry ("Peierls fitting exponents p", "0.5",
                           Patterns::List(Patterns::Double(0)),
                           "List of the first fitting exponents, $p$, between stress $\\sigma$ and the Peierls "
                           "stress $\\sigma_{\\text{peierls}}$ for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value. Units: none");
        prm.declare_entry ("Peierls fitting exponents q", "0.5",
                           Patterns::List(Patterns::Double(0)),
                           "List of the second fitting exponents, $q$, between stress $\\sigma$ and the Peierls "
                           "stress $\\sigma_{\\text{peierls}}$ for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value. Units: none");
      }



      template <int dim>
      void
      PeierlsCreep<dim>::parse_parameters (ParameterHandler &prm)
      {
        // increment by one for background:
        const unsigned int n_fields = this->n_compositional_fields() + 1;

        if (prm.get ("Peierls creep flow law") == "viscosity approximation")
          peierls_creep_flow_law = viscosity_approx;
        else
          AssertThrow(false, ExcMessage("Not a valid Peierls creep flow law"));

        prefactors_peierls = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Prefactors for Peierls creep"))),
                                                                     n_fields,
                                                                     "Prefactors for Peierls creep");

        stress_exponents_peierls = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Stress exponents for Peierls creep"))),
                                                                           n_fields,
                                                                           "Stress exponents for Peierls creep");

        activation_energies_peierls = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Activation energies for Peierls creep"))),
                                                                              n_fields,
                                                                              "Activation energies for Peierls creep");

        activation_volumes_peierls = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Activation volumes for Peierls creep"))),
                                                                             n_fields,
                                                                             "Activation volumes for Peierls creep");

        peierls_stresses = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Peierls stresses"))),
                                                                   n_fields,
                                                                   "Peierls stresses");
        peierls_fitting_parameters = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Peierls fitting parameters"))),
                                                                             n_fields,
                                                                             "Peierls fitting parameters");

        peierls_fitting_exponents_p = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Peierls fitting exponents p"))),
                                                                              n_fields,
                                                                              "Peierls fitting exponents p");

        peierls_fitting_exponents_q = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Peierls fitting exponents q"))),
                                                                              n_fields,
                                                                              "Peierls fitting exponents q");
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
  }
}
