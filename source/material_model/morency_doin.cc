/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <aspect/material_model/morency_doin.h>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    std::vector<double>
    MorencyDoin<dim>::
    compute_volume_fractions( const std::vector<double> &compositional_fields) const
    {
      std::vector<double> volume_fractions( compositional_fields.size()+1);

      //clip the compositional fields so they are between zero and one
      std::vector<double> x_comp = compositional_fields;
      for ( unsigned int i=0; i < x_comp.size(); ++i)
        x_comp[i] = std::min(std::max(x_comp[i], 0.0), 1.0);

      //sum the compositional fields for normalization purposes
      double sum_composition = 0.0;
      for ( unsigned int i=0; i < x_comp.size(); ++i)
        sum_composition += x_comp[i];

      if (sum_composition >= 1.0)
        {
          volume_fractions[0] = 0.0;  //background mantle
          for ( unsigned int i=1; i <= x_comp.size(); ++i)
            volume_fractions[i] = x_comp[i-1]/sum_composition;
        }
      else
        {
          volume_fractions[0] = 1.0 - sum_composition; //background mantle
          for ( unsigned int i=1; i <= x_comp.size(); ++i)
            volume_fractions[i] = x_comp[i-1];
        }
      return volume_fractions;
    }

    template <int dim>
    void
    MorencyDoin<dim>::
    evaluate(const MaterialModelInputs &in,
             MaterialModelOutputs &out) const
    {
      const double R = 8.32; // J mol-1 K-1
      const double g = 9.8; // m s-2
      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          const Point<dim> position = in.position[i];
          const double temperature = in.temperature[i];
          const double pressure = in.pressure[i];
          const std::vector<double> composition = in.composition[i];
          const std::vector<double> volume_fractions = compute_volume_fractions(composition);
          const SymmetricTensor<2,dim> strain_rate = in.strain_rate[i];

          const double z = this->get_geometry_model().depth(position); // units: m

          double activation_energy = 0.0;
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            activation_energy += volume_fractions[j] * activation_energies[j] * 1e3; // Converted to J/mol to make units work
          double nv = 0.0;
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            nv += volume_fractions[j] * nvs[j];
          double np = 0.0;
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            np += volume_fractions[j] * nps[j];

          const double tauy = tau_0 + gamma*reference_density()*g*z;
          const double e2inv = std::sqrt(std::pow(second_invariant(strain_rate),2) + std::pow(min_strain_rate,2));
          const double veffv = B * std::pow(e2inv/ref_strain_rate, -1.0+1.0/nv) * std::exp((activation_energy+activation_volume*reference_density()*g*z)/(nv*R*temperature));
          const double veffp = tauy * (std::pow(e2inv, -1.0+1.0/np) / std::pow(ref_strain_rate, 1.0/np));
          // Effective viscosity = harmonic mean of diffusion and dislocation creep. Range is limited to 1e17-1e28 for stability.
          const double veff = std::min(std::max(veff_coefficient * std::pow((1.0/veffp + 1.0/veffv), -1.0), min_visc), max_visc);
          out.viscosities[i] = veff;

          double density = 0.0;
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            {
              //not strictly correct if thermal expansivities are different, since we are interpreting
              //these compositions as volume fractions, but the error introduced should not be too bad.
              const double temperature_factor= (1.0 - thermal_expansivities[j] * (temperature - reference_T));
              density += volume_fractions[j] * densities[j] * temperature_factor;
            }
          out.densities[i] = density;

          // Thermal expansion coefficients at the given positions.
          double thermal_expansivity = 0.0;
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            thermal_expansivity += volume_fractions[j] * thermal_expansivities[j];
          out.thermal_expansion_coefficients[i] = thermal_expansivity;
          // Specific heat at the given positions.
          out.specific_heat[i] = heat_capacity;
          // Thermal conductivity at the given positions.
          out.thermal_conductivities[i] = thermal_diffusivity * heat_capacity * density;
          // Compressibility at the given positions.
          // The compressibility is given as
          // $\frac 1\rho \frac{\partial\rho}{\partial p}$.
          out.compressibilities[i] = 0.0;
          // Pressure derivative of entropy at the given positions.
          out.entropy_derivative_pressure[i] = 0.0;
          // Temperature derivative of entropy at the given positions.
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
    MorencyDoin<dim>::
    reference_viscosity () const
    {
      return ref_visc;
    }

    template <int dim>
    double
    MorencyDoin<dim>::
    reference_density () const
    {
      return densities[0];
    }

    template <int dim>
    bool
    MorencyDoin<dim>::
    viscosity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      // compare this with the implementation of the viscosity() function
      // to see the dependencies
      if (((dependence & NonlinearDependence::temperature) != NonlinearDependence::none))
        return true;
      else if ((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none)
        return true;
      else if (((dependence & NonlinearDependence::pressure) != NonlinearDependence::none))
        return true;
      else if (((dependence & NonlinearDependence::strain_rate) != NonlinearDependence::none))
        return true;
      else
        return false;
    }

    template <int dim>
    bool
    MorencyDoin<dim>::
    density_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      // compare this with the implementation of the density() function
      // to see the dependencies
      if ((dependence & NonlinearDependence::temperature) != NonlinearDependence::none)
        return true;
      else if ((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none)
        return true;
      else if (((dependence & NonlinearDependence::pressure) != NonlinearDependence::none))
        return true;
      else
        return false;
    }

    template <int dim>
    bool
    MorencyDoin<dim>::
    compressibility_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    MorencyDoin<dim>::
    specific_heat_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    MorencyDoin<dim>::
    thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      if ((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none)
        return true;
      else
        return false;
    }

    template <int dim>
    bool
    MorencyDoin<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    MorencyDoin<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Morency and Doin");
        {
          prm.declare_entry ("Densities", "3300.",
                             Patterns::List(Patterns::Double(0)),
                             "List of densities, $\\rho$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: $kg / m^3$");
          prm.declare_entry ("Activation energies", "500",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation energies, $E_a$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: $kJ / mol$");
          prm.declare_entry ("Thermal expansivities", "3.5e-5",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal expansivities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: $1 / K$");
          prm.declare_entry ("Stress exponents for viscous rheology", "3",
                             Patterns::List(Patterns::Double(0)),
                             "List of stress exponents, $n_v$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: None");
          prm.declare_entry ("Stress exponents for plastic rheology", "30",
                             Patterns::List(Patterns::Double(0)),
                             "List of stress exponents, $n_p$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: None");
          prm.declare_entry ("Thermal diffusivity", "0.8e-6", Patterns::Double(0), "Units: $m^2/s$");
          prm.declare_entry ("Heat capacity", "1.25e3", Patterns::Double(0), "Units: $J / (K * kg)$");
          prm.declare_entry ("Activation volume", "6.4e-6", Patterns::Double(0), "($V_a$). Units: $m^3 / mol$");
          prm.declare_entry ("Reference strain rate", "6.4e-16", Patterns::Double(0), "($\\dot{\\varepsilon}_{ref}$). Units: $1 / s$");
          prm.declare_entry ("Preexponential constant for viscous rheology law", "1.24e14", Patterns::Double(0), "($B$). Units: None");
          prm.declare_entry ("Coefficient of yield stress increase with depth", "0.25", Patterns::Double(0), "($\\gamma$). Units: None");
          prm.declare_entry ("Cohesive strength of rocks at the surface", "117", Patterns::Double(0), "($\\tau_0$). Units: $Pa$");
          prm.declare_entry ("Reference temperature", "293", Patterns::Double(0), "For calculating density by thermal expansivity. Units: $K$");
          prm.declare_entry ("Minimum strain rate", "1.4e-20", Patterns::Double(0), "Stabilizes strain dependent viscosity. Units: $1 / s$");

          prm.declare_entry ("Minimum viscosity", "1e17", Patterns::Double(0), "Lower cutoff for effective viscosity. Units: $Pa s$");
          prm.declare_entry ("Maximum viscosity", "1e28", Patterns::Double(0), "Upper cutoff for effective viscosity. Units: $Pa s$");
          prm.declare_entry ("Effective viscosity coefficient", "1.0", Patterns::Double(0), "Scaling coefficient for effective viscosity.");
          prm.declare_entry ("Reference viscosity", "1e22", Patterns::Double(0), "Reference viscosity for nondimensionalization.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    MorencyDoin<dim>::parse_parameters (ParameterHandler &prm)
    {
      //increment by one for background:
      const unsigned int n_fields = this->n_compositional_fields() + 1;

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Morency and Doin");
        {
          gamma = prm.get_double("Coefficient of yield stress increase with depth");
          thermal_diffusivity = prm.get_double("Thermal diffusivity");
          heat_capacity = prm.get_double("Heat capacity");
          activation_volume = prm.get_double("Activation volume");
          ref_strain_rate = prm.get_double("Reference strain rate");
          B = prm.get_double("Preexponential constant for viscous rheology law");
          tau_0 = prm.get_double ("Cohesive strength of rocks at the surface");
          min_strain_rate = prm.get_double("Minimum strain rate");
          reference_T = prm.get_double("Reference temperature");

          min_visc = prm.get_double ("Minimum viscosity");
          max_visc = prm.get_double ("Maximum viscosity");
          veff_coefficient = prm.get_double ("Effective viscosity coefficient");
          ref_visc = prm.get_double ("Reference viscosity");

          std::vector<double> x_values;

          // Parse densities
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Densities")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of density list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            densities.assign( n_fields , x_values[0]);
          else
            densities = x_values;

          // Parse activation energies
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Activation energies")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of activation energy list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            activation_energies.assign( n_fields , x_values[0] );
          else
            activation_energies = x_values;

          // Parse thermal expansivities
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Thermal expansivities")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of thermal expansivity list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            thermal_expansivities.assign( n_fields , x_values[0]);
          else
            thermal_expansivities = x_values;

          // Parse stress exponents
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Stress exponents for viscous rheology")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of nv list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            nvs.assign( n_fields , x_values[0]);
          else
            nvs = x_values;
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Stress exponents for plastic rheology")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of np list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            nps.assign( n_fields , x_values[0]);
          else
            nps = x_values;
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    ASPECT_REGISTER_MATERIAL_MODEL(MorencyDoin,
                                   "Morency and Doin",
                                   "An implementation of the visco-plastic rheology described by (Morency"
                                   " and Doin, 2004). Compositional fields can each be assigned individual"
                                   " activation energies, reference densities, thermal expansivities,"
                                   " and stress exponents. The effective viscosity is defined as"
                                   "\n\n"
                                   " \\[v_{eff} = \\left(\\frac{1}{v_{eff}^v}+\\frac{1}{v_{eff}^p}\\right)^{-1}\\]"
                                   " where"
                                   " \\[v_{eff}^v = B \\left(\\frac{\\dot{\\varepsilon}}{\\dot{\\varepsilon}_{ref}}\\right)^{-1+1/n_v}"
                                   " exp\\left(\\frac{E_a +V_a \\rho_m g z}{n_v R T}\\right) \\]"
                                   " \\[v_{eff}^p = (\\tau_0 + \\gamma \\rho_m g z) \\left( \\frac{\\dot{\\varepsilon}^{-1+1/n_p}}"
                                   " {\\dot{\\varepsilon}_{ref}^{1/n_p}} \\right) \\]"
                                   "\n\n"
                                   " Where $B$ is a scaling constant, $\\dot{\\varepsilon}$ is"
                                   " related to the second invariant of the strain rate tensor, $\\dot{\\varepsilon}_{ref}$ is a"
                                   " reference strain rate, $n_v$ and $n_p$ are stress exponents, $E_a$ is the activation energy,"
                                   " $V_a$ is the activation volume, $\\rho_m$ is the mantle density, $R$ is the gas constant, $T$"
                                   " is temperature, $\\tau_0$ is the cohestive strength of rocks at the surface, $\\gamma$ is a"
                                   " coefficient of yield stress increase with depth, and $z$ is depth."
                                   " \n\n"
                                   " Morency, C., and M‐P. Doin. \"Numerical simulations of the mantle lithosphere delamination.\""
                                   " Journal of Geophysical Research: Solid Earth (1978–2012) 109.B3 (2004)."
                                   "\n\n"
                                   " The value for the components of this formula and additional"
                                   " parameters are read from the parameter file in subsection"
                                   " 'Material model/Morency and Doin'.")
  }
}
