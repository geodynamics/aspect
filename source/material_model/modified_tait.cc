/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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
    initialize()
    {
      const double Kdprime_0 = -reference_Kprime/reference_isothermal_bulk_modulus;

      a = (1. + reference_Kprime) / (1. + reference_Kprime  + reference_isothermal_bulk_modulus * Kdprime_0);
      b = reference_Kprime  / reference_isothermal_bulk_modulus - Kdprime_0  / (1. + reference_Kprime );
      c = (1. + reference_Kprime  + reference_isothermal_bulk_modulus * Kdprime_0) /
          (reference_Kprime * reference_Kprime + reference_Kprime - reference_isothermal_bulk_modulus * Kdprime_0);
    }

    template <int dim>
    void
    ModifiedTait<dim>::
    evaluate(const MaterialModelInputs<dim> &in,
             MaterialModelOutputs<dim> &out) const
    {
      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          const double temperature = std::max(1e-10,in.temperature[i]);
          const double pressure = in.pressure[i];

          out.viscosities[i] = eta;
          out.thermal_conductivities[i] = k_value;

          // Einstein equations for the heat capacity at constant volume
          const double x_einstein = einstein_temperature/temperature;
          const double x_einstein0 = einstein_temperature/reference_temperature;

          // The following Einstein energies and heat capacities are divided through by 3nR.
          // This doesn't matter, as the equation of state relies only on ratios of these quantities
          double E_th0, E_th, C_V0, C_V;
          if (x_einstein0 < 35)
            {
              E_th0 = einstein_temperature * (0.5 + 1. / (std::exp(x_einstein0) - 1.0));
              C_V0 = x_einstein0 * x_einstein0 * std::exp(x_einstein0) / std::pow(std::exp(x_einstein0) - 1.0, 2.0);
            }
          else
            {
              E_th0 = einstein_temperature * 0.5;
              C_V0 = x_einstein0 * x_einstein0 * std::exp(-x_einstein0);
            }

          if (x_einstein < 35)
            {
              E_th = einstein_temperature * (0.5 + 1. / (std::exp(x_einstein) - 1.0));
              C_V = x_einstein * x_einstein * std::exp(x_einstein) / std::pow(std::exp(x_einstein) - 1.0, 2.0);
            }
          else
            {
              E_th = einstein_temperature * 0.5;
              C_V = x_einstein * x_einstein * std::exp(-x_einstein);
            }

          // The relative thermal pressure
          const double Pth = reference_thermal_expansivity * reference_isothermal_bulk_modulus * E_th / C_V0;
          const double Pth_rel = reference_thermal_expansivity * reference_isothermal_bulk_modulus * (E_th - E_th0) / C_V0;
          const double psubpth = pressure - reference_pressure - Pth_rel;

          // x = rho0/rho
          const double x = 1 - a * (1. - std::pow((1. + b * psubpth), -1.0 * c));

          // xi
          double xi = (C_V / C_V0);

          // Here we calculate the pressure effect on the heat capacity. It is a bit involved.
          const double dintVdpdT = (reference_thermal_expansivity * reference_isothermal_bulk_modulus / reference_rho * a * xi) * (
                                     std::pow((1. + b * psubpth), - c) - std::pow((1. - b * Pth), - c));

          const double dSdT = reference_isothermal_bulk_modulus / reference_rho * std::pow((xi * reference_thermal_expansivity), 2) * \
                              (std::pow((1. + b * psubpth), -1. - c) - std::pow((1. - b * Pth), -1. - c)) + \
                              dintVdpdT * (( 1 - 2./x + 2./(std::exp(x) - 1.) ) * einstein_temperature/(temperature*temperature));


          const Point<1> Tpoint(temperature);
          double rho = reference_rho / x;
          double isothermal_compressibility = 1./(reference_isothermal_bulk_modulus * (1. + b * psubpth) * (a + (1. - a) * std::pow((1. + b * psubpth), c)));
          double alpha = reference_thermal_expansivity * xi * 1. / ((1. + b * psubpth) * (a + (1. - a) * std::pow((1 + b * psubpth), c)));
          double heat_capacity = reference_heat_capacity_function.value(Tpoint) + temperature * dSdT;

          // The following line might be useful for some compressibility models
          // double isentropic_compressibility = isothermal_compressibility - alpha*alpha*temperature/(rho*heat_capacity);
          out.densities[i] = rho;
          out.compressibilities[i] = isothermal_compressibility;
          out.thermal_expansion_coefficients[i] = alpha;
          out.specific_heat[i] = heat_capacity;

          out.entropy_derivative_pressure[i] = 0.0;
          out.entropy_derivative_temperature[i] = 0.0;
          // Change in composition due to chemical reactions at the
          // given positions. The term reaction_terms[i][c] is the
          // change in compositional field c at point i.
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;
        }
    }

    template <int dim>
    double
    ModifiedTait<dim>::
    reference_viscosity () const
    {
      return eta;
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
                             Patterns::Double (0),
                             "Reference pressure $\\P_0$. Units: $Pa$.");
          prm.declare_entry ("Reference temperature", "298.15",
                             Patterns::Double (0),
                             "Reference temperature $\\T_0$. Units: $K$.");
          prm.declare_entry ("Reference density", "3300",
                             Patterns::Double (0),
                             "The density at the reference pressure and temperature. "
                             "Units: $kg/m^3$.");
          prm.declare_entry ("Reference isothermal bulk modulus", "125e9",
                             Patterns::Double (0),
                             "The isothermal bulk modulus at the reference pressure and temperature. "
                             "Units: $kg/m^3$.");
          prm.declare_entry ("Reference bulk modulus derivative", "4",
                             Patterns::Double (0),
                             "The value of the first pressure derivative of the isothermal bulk modulus "
                             "at the reference pressure and temperature. "
                             "Units: $kg/m^3$.");
          prm.declare_entry ("Reference thermal expansivity", "2e-5",
                             Patterns::Double (0),
                             "The thermal expansion coefficient at the reference pressure and temperature. "
                             "Units: $1/K$.");
          prm.declare_entry ("Einstein temperature", "600",
                             Patterns::Double (0),
                             "The Einstein temperature at the reference pressure and temperature. "
                             "Units: $K$.");
          prm.declare_entry ("Viscosity", "1e21",
                             Patterns::Double (0),
                             "The value of the constant viscosity $\\eta_0$. Units: $kg/m/s$.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0),
                             "The value of the constant thermal conductivity $k$. "
                             "Units: $W/m/K$.");

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
                                   "equation of state as written in the paper of "
                                   "Holland and Powell, 2011 "
                                   "(``An improved and extended internally consistent "
                                   "thermodynamic dataset for phases of petrological interest, "
                                   "involving a new equation of state for solids'', "
                                   "J. metamorphic Geol., 2011, 29, 333-383, "
                                   "\\url{http://dx.doi.org/10.1111/j.1525-1314.2010.00923.x}). "
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
