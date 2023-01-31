/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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

#include <aspect/material_model/modified_tait.h>
#include <aspect/adiabatic_conditions/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    ModifiedTait<dim>::
    evaluate(const MaterialModelInputs<dim> &in,
             MaterialModelOutputs<dim> &out) const
    {
      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          const double temperature = in.temperature[i];
          const double pressure = in.pressure[i];

          AssertThrow(temperature > 100., ExcMessage("The current temperature is <100 K. "
                                                     "The modified Tait equations of state is "
                                                     "not designed for such temperatures. "
                                                     "If this exception is triggered at the "
                                                     "beginning of a run, you may need to "
                                                     "set Adiabatic surface temperature "
                                                     "to a value > 100 in the parameter file."));
          out.viscosities[i] = eta;
          out.thermal_conductivities[i] = k_value;

          // Einstein equations for the heat capacity at constant volume
          const double x_einstein = einstein_temperature/temperature;
          const double x_einstein0 = einstein_temperature/reference_temperature;


          AssertThrow(x_einstein < 35., ExcMessage("The current temperature is less than 1/35th "
                                                   "of the Einstein temperature. This is not "
                                                   "expected for geological materials."));
          AssertThrow(x_einstein0 < 35., ExcMessage("The reference temperature of the "
                                                    "equation of state is less than 1/35th "
                                                    "of the Einstein temperature. This is not "
                                                    "expected for geological materials."));

          // The following Einstein energies and heat capacities are divided through by 3nR.
          // This doesn't matter, as the equation of state relies only on ratios of these quantities.
          const double E_th0 = einstein_temperature * (0.5 + 1. / (std::exp(x_einstein0) - 1.0));
          const double C_V0 = x_einstein0 * x_einstein0 * std::exp(x_einstein0) / std::pow(std::exp(x_einstein0) - 1.0, 2.0);

          const double E_th = einstein_temperature * (0.5 + 1. / (std::exp(x_einstein) - 1.0));
          const double C_V = x_einstein * x_einstein * std::exp(x_einstein) / std::pow(std::exp(x_einstein) - 1.0, 2.0);

          // The relative thermal pressure
          const double Pth_rel = reference_thermal_expansivity * reference_isothermal_bulk_modulus * (E_th - E_th0) / C_V0;
          const double psubpth = pressure - reference_pressure - Pth_rel;

          // x = rho0/rho
          const double x = 1 - tait_a * (1. - std::pow((1. + tait_b * psubpth), -1.0 * tait_c));

          // xi
          const double xi = (C_V / C_V0);

          // Here we calculate the pressure effect on the heat capacity. It is a bit involved.
          const double dintVdpdT = ((reference_thermal_expansivity *
                                     reference_isothermal_bulk_modulus / reference_rho *
                                     tait_a * xi) *
                                    (std::pow((1. + tait_b * psubpth), - tait_c) -
                                     std::pow((1. - tait_b * Pth_rel), - tait_c)));

          const double dSdT = (reference_isothermal_bulk_modulus / reference_rho *
                               std::pow((xi * reference_thermal_expansivity), 2) *
                               (std::pow((1. + tait_b * psubpth), -1. - tait_c) -
                                std::pow((1. - tait_b * Pth_rel), -1. - tait_c)) +
                               dintVdpdT * (( 1 - 2./x_einstein + 2./(std::exp(x_einstein) - 1.))
                                            * einstein_temperature/(temperature*temperature)));

          const double rho = reference_rho / x;
          const double isothermal_compressibility = 1./(reference_isothermal_bulk_modulus * (1. + tait_b * psubpth) * (tait_a + (1. - tait_a) * std::pow((1. + tait_b * psubpth), tait_c)));
          const double alpha = reference_thermal_expansivity * xi * 1. / ((1. + tait_b * psubpth) * (tait_a + (1. - tait_a) * std::pow((1 + tait_b * psubpth), tait_c)));
          const double heat_capacity = reference_heat_capacity_function.value(Point<1>(temperature)) + temperature * dSdT;

          // The following line might be useful for some compressibility models
          // double isentropic_compressibility = isothermal_compressibility - alpha*alpha*temperature/(rho*heat_capacity);
          out.densities[i] = rho;
          out.compressibilities[i] = isothermal_compressibility;
          out.thermal_expansion_coefficients[i] = alpha;
          out.specific_heat[i] = heat_capacity;

          out.entropy_derivative_pressure[i] = -alpha / rho;
          out.entropy_derivative_temperature[i] = heat_capacity/temperature;
          // Change in composition due to chemical reactions at the
          // given positions. The term reaction_terms[i][c] is the
          // change in compositional field c at point i.
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;
        }
    }



    template <int dim>
    bool
    ModifiedTait<dim>::
    is_compressible () const
    {
      return true;
    }



    template <int dim>
    void
    ModifiedTait<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Modified Tait model");
        {
          prm.declare_entry ("Reference pressure", "1e5",
                             Patterns::Double (0.),
                             "Reference pressure $P_0$. Units: \\si{\\pascal}.");
          prm.declare_entry ("Reference temperature", "298.15",
                             Patterns::Double (0.),
                             "Reference temperature $T_0$. Units: \\si{\\kelvin}.");
          prm.declare_entry ("Reference density", "3300.",
                             Patterns::Double (0.),
                             "The density at the reference pressure and temperature. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
          prm.declare_entry ("Reference isothermal bulk modulus", "125e9",
                             Patterns::Double (0.),
                             "The isothermal bulk modulus at the reference pressure and temperature. "
                             "Units: \\si{\\pascal}.");
          prm.declare_entry ("Reference bulk modulus derivative", "4.",
                             Patterns::Double (0.),
                             "The value of the first pressure derivative of the isothermal bulk modulus "
                             "at the reference pressure and temperature. "
                             "Units: None.");
          prm.declare_entry ("Reference thermal expansivity", "2e-5",
                             Patterns::Double (0.),
                             "The thermal expansion coefficient at the reference pressure and temperature. "
                             "Units: \\si{\\per\\kelvin}.");
          prm.declare_entry ("Einstein temperature", "600.",
                             Patterns::Double (0.),
                             "The Einstein temperature at the reference pressure and temperature. "
                             "Units: \\si{\\kelvin}.");
          prm.declare_entry ("Viscosity", "1e21",
                             Patterns::Double (0.),
                             "The value of the constant viscosity $\\eta_0$. "
                             "Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0.),
                             "The value of the constant thermal conductivity $k$. "
                             "Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");

          prm.enter_subsection("Reference heat capacity function");
          {
            Functions::ParsedFunction<1>::declare_parameters(prm,1);
            prm.declare_entry("Function expression","1.25e3");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    ModifiedTait<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Modified Tait model");
        {
          reference_pressure                = prm.get_double ("Reference pressure");
          reference_temperature             = prm.get_double ("Reference temperature");
          reference_rho                     = prm.get_double ("Reference density");
          reference_isothermal_bulk_modulus = prm.get_double ("Reference isothermal bulk modulus");
          reference_Kprime                  = prm.get_double ("Reference bulk modulus derivative");
          reference_thermal_expansivity     = prm.get_double ("Reference thermal expansivity");
          einstein_temperature              = prm.get_double ("Einstein temperature");

          eta                               = prm.get_double ("Viscosity");
          k_value                           = prm.get_double ("Thermal conductivity");

          prm.enter_subsection("Reference heat capacity function");
          {
            try
              {
                reference_heat_capacity_function.parse_parameters(prm);
              }
            catch (...)
              {
                std::cerr << "FunctionParser failed to parse\n"
                          << "\t Reference heat capacity function\n"
                          << "with expression \n"
                          << "\t' " << prm.get("Function expression") << "'";
                throw;
              }
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Calculate dependent parameters
      Kdprime_0 = -reference_Kprime/reference_isothermal_bulk_modulus;
      tait_a = (1. + reference_Kprime) / (1. + reference_Kprime  + reference_isothermal_bulk_modulus * Kdprime_0);
      tait_b = (reference_Kprime  / reference_isothermal_bulk_modulus) - (Kdprime_0  / (1. + reference_Kprime ));
      tait_c = (1. + reference_Kprime  + reference_isothermal_bulk_modulus * Kdprime_0) /
               (reference_Kprime * reference_Kprime + reference_Kprime - reference_isothermal_bulk_modulus * Kdprime_0);

      // Declare dependencies on solution variables
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;
      this->model_dependence.viscosity = NonlinearDependence::none;


      this->model_dependence.density = NonlinearDependence::temperature
                                       | NonlinearDependence::pressure
                                       | NonlinearDependence::compositional_fields;
      this->model_dependence.compressibility = NonlinearDependence::temperature
                                               | NonlinearDependence::pressure
                                               | NonlinearDependence::compositional_fields;
      this->model_dependence.specific_heat = NonlinearDependence::temperature
                                             | NonlinearDependence::pressure
                                             | NonlinearDependence::compositional_fields;

    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(ModifiedTait,
                                   "modified tait",
                                   "A material model that implements the thermal modified Tait "
                                   "equation of state as written in \\cite{HP2011}. "
                                   "Constant values are used for the thermal conductivity and "
                                   "viscosity. The defaults for all coefficients are chosen "
                                   "to be similar to what is believed to be correct "
                                   "for Earth's mantle. "
                                   "All of the values that define this model are read "
                                   "from a section ``Material model/Modified Tait model'' "
                                   "in the input file, see "
                                   "Section~\\ref{parameters:Material_20model/Modified_20Tait_20model}.")
  }
}
