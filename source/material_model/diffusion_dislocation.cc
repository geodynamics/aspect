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

#include <aspect/material_model/diffusion_dislocation.h>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    std::vector<double>
    DiffusionDislocation<dim>::
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
    double
    DiffusionDislocation<dim>::
    average_value ( const std::vector<double> &composition,
                    const std::vector<double> &parameter_values,
                    const enum averaging_scheme &average_type) const
    {
      double averaged_parameter = 0.0;
      const std::vector<double> volume_fractions = compute_volume_fractions(composition);

      switch (average_type)
        {
          case arithmetic:
          {
            for (unsigned int i=0; i< volume_fractions.size(); ++i)
              averaged_parameter += volume_fractions[i]*parameter_values[i];
            break;
          }
          case harmonic:
          {
            for (unsigned int i=0; i< volume_fractions.size(); ++i)
              averaged_parameter += volume_fractions[i]/(parameter_values[i]);
            averaged_parameter = 1.0/averaged_parameter;
            break;
          }
          case geometric:
          {
            for (unsigned int i=0; i < volume_fractions.size(); ++i)
              averaged_parameter += volume_fractions[i]*std::log(parameter_values[i]);
            averaged_parameter = std::exp(averaged_parameter);
            break;
          }
          case maximum_composition:
          {
            const unsigned int i = (unsigned int)(std::max_element( volume_fractions.begin(),
                                                                    volume_fractions.end() )
                                                  - volume_fractions.begin());
            averaged_parameter = parameter_values[i];
            break;
          }
          default:
          {
            AssertThrow( false, ExcNotImplemented() );
            break;
          }
        }
      return averaged_parameter;
    }


    template <int dim>
    std::vector<double>
    DiffusionDislocation<dim>::
    calculate_viscosities ( const std::vector<double> &volume_fractions,
                            const double &pressure,
                            const double &temperature,
                            const SymmetricTensor<2,dim> &strain_rate) const
    {
      // Viscosities
      const double e2inv = second_invariant(strain_rate);

      // ---- Find effective viscosities for each of the individual phases
      std::vector<double> composition_viscosities(volume_fractions.size()); // viscosities should have same number of entries as compositional fields
      for (unsigned int j=0; j < volume_fractions.size(); ++j)
        {
          const double viscosity_diffusion = std::min(1e22,(1e0/prefactors_diffusion[j])*
                                                      std::exp((activation_energies_diffusion[j]+activation_volumes_diffusion[j]*pressure)
                                                               /(constants::gas_constant*temperature)));

          double one_over_viscosity_dislocation = 0.0;
          if (e2inv > 2.0*std::numeric_limits<double>::min())
            one_over_viscosity_dislocation = (constants::gas_constant*temperature)
                                             /
                                             std::min(1e22,std::pow(prefactors_dislocation[j],-1e0/stress_exponents_dislocation[j])*
                                                      std::pow(e2inv,(1e0-stress_exponents_dislocation[j])/
                                                               stress_exponents_dislocation[j])*
                                                      std::exp((activation_energies_dislocation[j]+
                                                                activation_volumes_dislocation[j]*pressure)/(stress_exponents_dislocation[j])));

          composition_viscosities[j] = std::min(std::max(std::pow((1.0/viscosity_diffusion + one_over_viscosity_dislocation), -1.0), min_visc), max_visc);
        }
      return composition_viscosities;
    }

    template <int dim>
    void
    DiffusionDislocation<dim>::
    evaluate(const MaterialModelInputs &in,
             MaterialModelOutputs &out) const
    {
      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          // const Point<dim> position = in.position[i];
          const double temperature = in.temperature[i];
          const double pressure= in.pressure[i];
          const std::vector<double> composition = in.composition[i];
          const std::vector<double> volume_fractions = compute_volume_fractions(composition);

          // Averaging composition-field dependent properties

          // densities
          double density = 0.0;
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            {
              //not strictly correct if thermal expansivities are different, since we are interpreting
              //these compositions as volume fractions, but the error introduced should not be too bad.
              const double temperature_factor= (1.0 - thermal_expansivities[j] * (temperature - reference_T));
              density += volume_fractions[j] * densities[j] * temperature_factor;
            }

          // thermal expansivities
          double thermal_expansivity = 0.0;
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            thermal_expansivity += volume_fractions[j] * thermal_expansivities[j];

          // calculate effective viscosity
          if (in.strain_rate.size())
            {
              const std::vector<double> composition_viscosities = calculate_viscosities(volume_fractions, pressure, temperature, in.strain_rate[i]);
              const double veff = average_value(composition, composition_viscosities, viscosity_averaging);
              out.viscosities[i] = veff;
            }

          out.densities[i] = density;
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
    DiffusionDislocation<dim>::
    reference_viscosity () const
    {
      return ref_visc;
    }

    template <int dim>
    double
    DiffusionDislocation<dim>::
    reference_density () const
    {
      return densities[0];
    }

    template <int dim>
    bool
    DiffusionDislocation<dim>::
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
    DiffusionDislocation<dim>::
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
    DiffusionDislocation<dim>::
    compressibility_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    DiffusionDislocation<dim>::
    specific_heat_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    DiffusionDislocation<dim>::
    thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      if ((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none)
        return true;
      else
        return false;
    }

    template <int dim>
    bool
    DiffusionDislocation<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    DiffusionDislocation<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("DiffusionDislocation");
        {
          // Reference and minimum/maximum values
          prm.declare_entry ("Reference temperature", "293", Patterns::Double(0),
                             "For calculating density by thermal expansivity. Units: $K$");
          prm.declare_entry ("Minimum strain rate", "1.4e-20", Patterns::Double(0),
                             "Stabilizes strain dependent viscosity. Units: $1 / s$");
          prm.declare_entry ("Minimum viscosity", "1e17", Patterns::Double(0),
                             "Lower cutoff for effective viscosity. Units: $Pa s$");
          prm.declare_entry ("Maximum viscosity", "1e28", Patterns::Double(0),
                             "Upper cutoff for effective viscosity. Units: $Pa s$");
          prm.declare_entry ("Effective viscosity coefficient", "1.0", Patterns::Double(0),
                             "Scaling coefficient for effective viscosity.");
          prm.declare_entry ("Reference viscosity", "1e22", Patterns::Double(0),
                             "Reference viscosity for nondimensionalization. Units $Pa s$");


          // Equation of state parameters
          prm.declare_entry ("Thermal diffusivity", "0.8e-6", Patterns::Double(0), "Units: $m^2/s$");
          prm.declare_entry ("Heat capacity", "1.25e3", Patterns::Double(0), "Units: $J / (K * kg)$");
          prm.declare_entry ("Densities", "3300.",
                             Patterns::List(Patterns::Double(0)),
                             "List of densities, $\\rho$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: $kg / m^3$");
          prm.declare_entry ("Thermal expansivities", "3.5e-5",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal expansivities for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: $1 / K$");

          // Rheological parameters
          prm.declare_entry ("Viscosity averaging scheme", "harmonic",
                             Patterns::Selection("arithmetic|harmonic|geometric|maximum composition"),
                             "When more than one compositional field is present at a point "
                             "with different viscosities, we need to come up with an average "
                             "viscosity at that point.  Select a weighted harmonic, arithmetic, "
                             "geometric, or maximum composition.");
          // Diffusion creep parameters
          prm.declare_entry ("Prefactors for diffusion creep", "1.92e-11",
                             Patterns::List(Patterns::Double(0)),
                             "List of viscosity prefactors, $A$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: $1 / s$");
          prm.declare_entry ("Activation energies for diffusion creep", "335e3",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation energies, $E_a$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: $J / mol$");
          prm.declare_entry ("Activation volumes for diffusion creep", "6.4e-6",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation volumes, $V_a$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value.  Units: $m^3 / mol$");

          // Dislocation creep parameters
          prm.declare_entry ("Prefactors for dislocation creep", "2.42e-10",
                             Patterns::List(Patterns::Double(0)),
                             "List of viscosity prefactors, $A$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: $Pa^{-n_dislocation} s^{-1}$");
          prm.declare_entry ("Activation energies for dislocation creep", "540e3",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation energies, $E_a$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: $J / mol$");
          prm.declare_entry ("Activation volumes for dislocation creep", "6.4e-6",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation volumes, $V_a$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one value is given, then all use the same value.  Units: $m^3 / mol$");
          prm.declare_entry ("Stress exponents for dislocation creep", "3",
                             Patterns::List(Patterns::Double(0)),
                             "List of stress exponents, $n_dislocation$, for background mantle and compositional fields,"
                             "for a total of N+1 values, where N is the number of compositional fields."
                             "If only one values is given, then all use the same value.  Units: None");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    DiffusionDislocation<dim>::parse_parameters (ParameterHandler &prm)
    {
      //increment by one for background:
      const unsigned int n_fields = this->n_compositional_fields() + 1;

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("DiffusionDislocation");
        {
          // Initialise empty vector for compositional field variables
          std::vector<double> x_values;

          // Reference and minimum/maximum values
          reference_T = prm.get_double("Reference temperature");
          min_strain_rate = prm.get_double("Minimum strain rate");
          min_visc = prm.get_double ("Minimum viscosity");
          max_visc = prm.get_double ("Maximum viscosity");
          veff_coefficient = prm.get_double ("Effective viscosity coefficient");
          ref_visc = prm.get_double ("Reference viscosity");

          // Equation of state parameters
          thermal_diffusivity = prm.get_double("Thermal diffusivity");
          heat_capacity = prm.get_double("Heat capacity");

          // ---- Densities
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Densities")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of density list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            densities.assign( n_fields , x_values[0]);
          else
            densities = x_values;

          // ---- Thermal expansivities
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Thermal expansivities")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of thermal expansivity list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            thermal_expansivities.assign( n_fields , x_values[0]);
          else
            thermal_expansivities = x_values;


          // Rheological parameters
          if (prm.get ("Viscosity averaging scheme") == "harmonic")
            viscosity_averaging = harmonic;
          else if (prm.get ("Viscosity averaging scheme") == "arithmetic")
            viscosity_averaging = arithmetic;
          else if (prm.get ("Viscosity averaging scheme") == "geometric")
            viscosity_averaging = geometric;
          else if (prm.get ("Viscosity averaging scheme") == "maximum composition")
            viscosity_averaging = maximum_composition;
          else
            AssertThrow(false, ExcMessage("Not a valid viscosity averaging scheme"));

          // Diffusion creep parameters
          // ---- diffusion creep prefactors
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Prefactors for diffusion creep")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of prefactors for diffusion list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            prefactors_diffusion.assign( n_fields , x_values[0] );
          else
            prefactors_diffusion = x_values;

          // ---- diffusion creep activation energies
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Activation energies for diffusion creep")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of activation energy for diffusion list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            activation_energies_diffusion.assign( n_fields , x_values[0] );
          else
            activation_energies_diffusion = x_values;

          // ---- diffusion creep activation volumes
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Activation volumes for diffusion creep")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of activation volume for diffusion list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            activation_volumes_diffusion.assign( n_fields , x_values[0] );
          else
            activation_volumes_diffusion = x_values;


          // Dislocation creep parameters
          // ---- dislocation creep prefactors
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Prefactors for dislocation creep")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of prefactors for dislocation list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            prefactors_dislocation.assign( n_fields , x_values[0] );
          else
            prefactors_dislocation = x_values;

          // ---- dislocation creep activation energies
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Activation energies for dislocation creep")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of activation energy for dislocation list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            activation_energies_dislocation.assign( n_fields , x_values[0] );
          else
            activation_energies_dislocation = x_values;

          // ---- dislocation creep activation volumes
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Activation volumes for dislocation creep")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of activation volume for dislocation list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            activation_volumes_dislocation.assign( n_fields , x_values[0] );
          else
            activation_volumes_dislocation = x_values;

          // ---- dislocation creep stress exponents
          x_values = Utilities::string_to_double(Utilities::split_string_list(prm.get ("Stress exponents for dislocation creep")));
          AssertThrow(x_values.size() == 1u || (x_values.size() == n_fields),
                      ExcMessage("Length of stress exponents for dislocation list must be either one, or n_compositional_fields+1"));
          if (x_values.size() == 1)
            stress_exponents_dislocation.assign( n_fields , x_values[0] );
          else
            stress_exponents_dislocation = x_values;
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    ASPECT_REGISTER_MATERIAL_MODEL(DiffusionDislocation,
                                   "diffusion dislocation",
                                   " An implementation of a viscous rheology including diffusion"
                                   " and dislocation creep."
                                   " Compositional fields can each be assigned individual"
                                   " activation energies, reference densities, thermal expansivities,"
                                   " and stress exponents. The effective viscosity is defined as"
                                   "\n\n"
                                   " \\[v_\\text{eff} = \\left(\\frac{1}{v_\\text{eff}^\\text{diff}}+\\frac{1}{v_\\text{eff}^\\text{dis}\\right)^{-1}\\]"
                                   " where"
                                   " \\[v_\\text{eff}^\\text{diff} = A_\\text{diff}^{-1} \\exp\\left(\frac{E_\\text{diff} + PV_\\text{diff}}{RT}\\right)\\]"
                                   " \\[v_\\text{eff}^\\text{dis} =  A_\\text{dis}^{\\frac{-1}{n_{dis}}} \\dot{\\varepsilon}^{\frac{1-n}{n}} "
                                   "                                 \\exp\\left(\frac{E_\\text{diff} + PV_\\text{diff}}{n_\\text{dis}RT}\\right)\\]"
                                   "\n\n"
                                   " where $\\dot{\\varepsilon}$ is the second invariant of the strain rate tensor,"
                                   " $A_i$ are prefactors where $i$ corresponds to diffusion or dislocation creep,"
                                   " $E_i$ are the activation energies, $V_i$ are the activation volumes,"
                                   " $n_dislocation$ are stress exponents for dislocation creep,"
                                   " $\\rho_m$ is the mantle density, $R$ is the gas constant,"
                                   " $T$ is temperature, and $P$ is pressure."
                                   " \n\n"
                                   " The value for the components of this formula and additional"
                                   " parameters are read from the parameter file in subsection"
                                   " 'Material model/DiffusionDislocation'.")
  }
}
