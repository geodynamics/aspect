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
                                            const double temperature,
                                            const unsigned int composition) const
      {
        double viscosity = 0.0;

        switch (peierls_creep_flow_law)
          {
            case mei_viscosity_approx:
            {
              /**
               * An approximation of the Peierls creep formulation from Mei et al. (2010), where stress is
               * replaced with strain rate (second invariant) when calculating viscosity. This substitution
               * require two additional fitting parameters (gamma,p). Details of this derivation can be found
               * in the following document: https://ucdavis.app.box.com/s/scaorcblr9u294836pgk4hyf60gv7psr.
               * The formulation for the viscosity is:
               *  ((0.5*gamma*sigma_p)/(A*((gamma*sigma_p)^2)^(1/(s+2)))) * exp((E/(R*T))*((1-gamma^p)/(s+2))) * edot_ii^(1/(s+2)-1)
               * where
               * s = (E/(R*T))*p*(gamma^p)
               * sigma_p is the Peierls stress
               * gamma is the Peierls creep fitting parameter
               * p is the Peierls creep fitting exponent
               * A is the Peierls prefactor term (1/s 1/Pa^2)
               * E is the Peierls activation energy (J/mol)
               * edot_ii is the second invariant of the deviatoric strain rate
               * T is temperature (K)
               * R is the gas constant
               */
              const double s = (activation_energies_peierls[composition] / (constants::gas_constant * temperature)) *
                               peierls_fitting_exponents[composition] * std::pow(peierls_fitting_parameters[composition],peierls_fitting_exponents[composition]);

              const double left_term = 0.5 * (peierls_fitting_parameters[composition] * peierls_stresses[composition]) /
                                       std::pow((prefactors_peierls[composition] * std::pow(peierls_fitting_parameters[composition] * peierls_stresses[composition],2)),( 1 / (s + 2)));

              const double middle_term = std::exp((activation_energies_peierls[composition] / (constants::gas_constant * temperature)) *
                                                  ((1. - std::pow(peierls_fitting_parameters[composition],peierls_fitting_exponents[composition])) / (s + 2)));

              const double right_term = std::pow(strain_rate, (1. / (s + 2)) - 1.);

              viscosity = left_term * middle_term * right_term;
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
        prm.declare_entry ("Peierls creep flow law", "mei viscosity approximation",
                           Patterns::Selection("mei viscosity approximation"),
                           "Select what type of Peierls creep flow law to use. Currently, the "
                           "only available option is an approximation to the Mei 2010 Peierls "
                           "creep flow law, which uses the strain rate invariant, rather than "
                           "stress, as input. Future options will allow formulations that use "
                           "stress as an input, instead of strain rate.");
        prm.declare_entry ("Prefactors for Peierls creep", "1.4e-7",
                           Patterns::List(Patterns::Double(0)),
                           "List of viscosity prefactors, $A$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value. Units: $1/s 1/Pa^{2}$");
        prm.declare_entry ("Activation energies for Peierls creep", "320e3",
                           Patterns::List(Patterns::Double(0)),
                           "List of activation energies, $E$, for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value.  Units: $J / mol$");
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
        prm.declare_entry ("Peierls fitting exponents", "0.5",
                           Patterns::List(Patterns::Double(0)),
                           "List of fitting exponents, $p$, between stress $\\sigma$ and the Peierls stres "
                           "$\\sigma_{\\text{peierls}}$ for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value. Units: none");
      }



      template <int dim>
      void
      PeierlsCreep<dim>::parse_parameters (ParameterHandler &prm)
      {
        // increment by one for background:
        const unsigned int n_fields = this->n_compositional_fields() + 1;

        if (prm.get ("Peierls creep flow law") == "mei viscosity approximation")
          peierls_creep_flow_law = mei_viscosity_approx;
        else
          AssertThrow(false, ExcMessage("Not a valid Peierls creep flow law"));

        prefactors_peierls = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Prefactors for Peierls creep"))),
                                                                     n_fields,
                                                                     "Prefactors for Peierls creep");
        activation_energies_peierls = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Activation energies for Peierls creep"))),
                                                                              n_fields,
                                                                              "Activation energies for Peierls creep");
        peierls_stresses = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Peierls stresses"))),
                                                                   n_fields,
                                                                   "Peierls stresses");
        peierls_fitting_parameters = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Peierls fitting parameters"))),
                                                                             n_fields,
                                                                             "Peierls fitting parameters");

        peierls_fitting_exponents = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Peierls fitting exponents"))),
                                                                            n_fields,
                                                                            "Peierls fitting exponents");
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
