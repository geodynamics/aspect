/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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

#include <aspect/material_model/visco_plastic.h>
#include <aspect/utilities.h>
using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    std::vector<double>
    ViscoPlastic<dim>::
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
          volume_fractions[0] = 0.0;  //background material
          for ( unsigned int i=1; i <= x_comp.size(); ++i)
            volume_fractions[i] = x_comp[i-1]/sum_composition;
        }
      else
        {
          volume_fractions[0] = 1.0 - sum_composition; //background material
          for ( unsigned int i=1; i <= x_comp.size(); ++i)
            volume_fractions[i] = x_comp[i-1];
        }
      return volume_fractions;
    }

    template <int dim>
    double
    ViscoPlastic<dim>::
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
    ViscoPlastic<dim>::
    calculate_isostrain_viscosities ( const std::vector<double> &volume_fractions,
                                      const double &pressure,
                                      const double &temperature,
                                      const SymmetricTensor<2,dim> &strain_rate,
                                      const ViscosityScheme &viscous_type,
                                      const YieldScheme &yield_type) const
    {
      // This function calculates viscosities assuming that all the compositional fields
      // experience the same strain rate (isostrain).

      // Calculate the square root of the second moment invariant for the deviatoric strain rate tensor.
      // The first time this function is called (first iteration of first time step)
      // a specified "reference" strain rate is used as the returned value would
      // otherwise be zero.
      const double edot_ii = ( (this->get_timestep_number() == 0 && strain_rate.norm() <= std::numeric_limits<double>::min())
                               ?
                               ref_strain_rate
                               :
                               std::max(std::sqrt(std::fabs(second_invariant(deviator(strain_rate)))),
                                        min_strain_rate) );

      // Calculate viscosities for each of the individual compositional phases
      std::vector<double> composition_viscosities(volume_fractions.size());
      for (unsigned int j=0; j < volume_fractions.size(); ++j)
        {
          // Power law creep equation
          //    viscosity = 0.5 * A^(-1/n) * edot_ii^((1-n)/n) * d^(-m/n) * exp((E + P*V)/(nRT))
          // A: prefactor, edot_ii: square root of second invariant of deviatoric strain rate tensor,
          // d: grain size, m: grain size exponent, E: activation energy, P: pressure,
          // V; activation volume, n: stress exponent, R: gas constant, T: temperature.
          // Note: values of A, d, m, E, V and n are distinct for diffusion & dislocation creep

          // Diffusion creep: viscosity is grain size dependent (m!=0) and strain-rate independent (n=1)
          double viscosity_diffusion = 0.5 * std::pow(prefactors_diffusion[j],-1/stress_exponents_diffusion[j]) *
                                       std::exp((activation_energies_diffusion[j] + pressure*activation_volumes_diffusion[j])/
                                                (constants::gas_constant*temperature*stress_exponents_diffusion[j])) *
                                       std::pow(grain_size, -grain_size_exponents_diffusion[j]);

          // For dislocation creep, viscosity is grain size independent (m=0) and strain-rate dependent (n>1)
          double viscosity_dislocation = 0.5 * std::pow(prefactors_dislocation[j],-1/stress_exponents_dislocation[j]) *
                                         std::exp((activation_energies_dislocation[j] + pressure*activation_volumes_dislocation[j])/
                                                  (constants::gas_constant*temperature*stress_exponents_dislocation[j])) *
                                         std::pow(edot_ii,((1. - stress_exponents_dislocation[j])/stress_exponents_dislocation[j]));

          // Composite viscosity
          double viscosity_composite = (viscosity_diffusion * viscosity_dislocation)/(viscosity_diffusion + viscosity_dislocation);

          // Select what form of viscosity to use (diffusion, dislocation or composite)
          double viscosity_pre_yield = 0.0;
          switch (viscous_type)
            {
              case diffusion:
              {
                viscosity_pre_yield = viscosity_diffusion;
                break;
              }
              case dislocation:
              {
                viscosity_pre_yield = viscosity_dislocation;
                break;
              }
              case composite:
              {
                viscosity_pre_yield = viscosity_composite;
                break;
              }
              default:
              {
                AssertThrow( false, ExcNotImplemented() );
                break;
              }
            }

          // Calculate viscous stress
          double viscous_stress = 2. * viscosity_pre_yield * edot_ii;

          // Convert friction angle from degrees to radians
          double phi = angles_internal_friction[j] * numbers::PI/180.0;

          // Calculate drucker prager yield strength (i.e. yield stress)
          double yield_strength = ( (dim==3)
                                    ?
                                    ( 6.0 * cohesions[j] * std::cos(phi) + 2.0 * std::max(pressure,0.0) * std::sin(phi) )
                                    / ( std::sqrt(3.0) * (3.0 + std::sin(phi) ) )
                                    :
                                    cohesions[j] * std::cos(phi) + std::max(pressure,0.0) * std::sin(phi) );



          // If the viscous stress is greater than the yield strength, rescale the viscosity back to yield surface
          double viscosity_drucker_prager;
          if ( viscous_stress >= yield_strength  )
            {
              viscosity_drucker_prager = yield_strength / (2.0 * edot_ii);
            }
          else
            {
              viscosity_drucker_prager = viscosity_pre_yield;
            }


          // Stress limiter rheology
          double viscosity_limiter;
          viscosity_limiter = yield_strength / (2.0 * ref_strain_rate) *
                              std::pow((edot_ii/ref_strain_rate), 1./exponents_stress_limiter[j] - 1.0);

          // Select if yield viscosity is based on drucker prager or stress limiter rheology
          double viscosity_yield;
          switch (yield_type)
            {
              case stress_limiter:
              {
                //viscosity_yield = std::min(viscosity_limiter, viscosity_pre_yield);
                viscosity_yield = viscosity_limiter;
                break;
              }
              case drucker_prager:
              {
                viscosity_yield = viscosity_drucker_prager;
                break;
              }
              default:
              {
                AssertThrow( false, ExcNotImplemented() );
                break;
              }
            }

          // Limit the viscosity with specified minimum and maximum bounds
          composition_viscosities[j] = std::min(std::max(viscosity_yield, min_visc), max_visc);


        }
      return composition_viscosities;
    }

    template <int dim>
    void
    ViscoPlastic<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          const double temperature = in.temperature[i];
          const double pressure = in.pressure[i];
          const std::vector<double> &composition = in.composition[i];
          const std::vector<double> volume_fractions = compute_volume_fractions(composition);

          // Averaging composition-field dependent properties

          // densities
          double density = 0.0;
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            {
              //not strictly correct if thermal expansivities are different, since we are interpreting
              //these compositions as volume fractions, but the error introduced should not be too bad.
              const double temperature_factor = (1.0 - thermal_expansivities[j] * (temperature - reference_T));
              density += volume_fractions[j] * densities[j] * temperature_factor;
            }

          // thermal expansivities
          double thermal_expansivity = 0.0;
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            thermal_expansivity += volume_fractions[j] * thermal_expansivities[j];

          // heat capacities
          double heat_capacity = 0.0;
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            heat_capacity += volume_fractions[j] * heat_capacities[j];

          // thermal diffusivities
          double thermal_diffusivity = 0.0;
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            thermal_diffusivity += volume_fractions[j] * thermal_diffusivities[j];

          // calculate effective viscosity
          if (in.strain_rate.size())
            {
              // Currently, the viscosities for each of the compositional fields are calculated assuming
              // isostrain amongst all compositions, allowing calculation of the viscosity ratio.
              // TODO: This is only consistent with viscosity averaging if the arithmetic averaging
              // scheme is chosen. It would be useful to have a function to calculate isostress viscosities.
              const std::vector<double> composition_viscosities =
                calculate_isostrain_viscosities(volume_fractions, pressure, temperature, in.strain_rate[i],viscous_flow_law,yield_mechanism);

              // The isostrain condition implies that the viscosity averaging should be arithmetic (see above).
              // We have given the user freedom to apply alternative bounds, because in diffusion-dominated
              // creep (where n_diff=1) viscosities are stress and strain-rate independent, so the calculation
              // of compositional field viscosities is consistent with any averaging scheme.
              out.viscosities[i] = average_value(composition, composition_viscosities, viscosity_averaging);

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
    ViscoPlastic<dim>::
    reference_viscosity () const
    {
      return ref_visc;
    }

    template <int dim>
    double
    ViscoPlastic<dim>::
    reference_density () const
    {
      return densities[0];
    }

    template <int dim>
    bool
    ViscoPlastic<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    ViscoPlastic<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Visco Plastic");
        {
          // Reference and minimum/maximum values
          prm.declare_entry ("Reference temperature", "293", Patterns::Double(0),
                             "For calculating density by thermal expansivity. Units: $K$");
          prm.declare_entry ("Minimum strain rate", "1.0e-20", Patterns::Double(0),
                             "Stabilizes strain dependent viscosity. Units: $1 / s$");
          prm.declare_entry ("Reference strain rate","1.0e-15",Patterns::Double(0),
                             "Reference strain rate for first time step. Units: $1 / s$");
          prm.declare_entry ("Minimum viscosity", "1e17", Patterns::Double(0),
                             "Lower cutoff for effective viscosity. Units: $Pa s$");
          prm.declare_entry ("Maximum viscosity", "1e28", Patterns::Double(0),
                             "Upper cutoff for effective viscosity. Units: $Pa s$");
          prm.declare_entry ("Reference viscosity", "1e22", Patterns::Double(0),
                             "Reference viscosity for nondimensionalization. Units $Pa s$");

          // Equation of state parameters
          prm.declare_entry ("Thermal diffusivities", "0.8e-6",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal diffusivities, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: $m^2/s$");
          prm.declare_entry ("Heat capacities", "1.25e3",
                             Patterns::List(Patterns::Double(0)),
                             "List of heat capacities $C_p$, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: $J/kg/K$");
          prm.declare_entry ("Densities", "3300.",
                             Patterns::List(Patterns::Double(0)),
                             "List of densities, $\\rho$, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: $kg / m^3$");
          prm.declare_entry ("Thermal expansivities", "3.5e-5",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal expansivities for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: $1 / K$");

          // Rheological parameters
          prm.declare_entry ("Grain size", "1e-3", Patterns::Double(0), "Units: $m$");
          prm.declare_entry ("Viscosity averaging scheme", "harmonic",
                             Patterns::Selection("arithmetic|harmonic|geometric|maximum composition"),
                             "When more than one compositional field is present at a point "
                             "with different viscosities, we need to come up with an average "
                             "viscosity at that point.  Select a weighted harmonic, arithmetic, "
                             "geometric, or maximum composition.");
          prm.declare_entry ("Viscous flow law", "composite",
                             Patterns::Selection("diffusion|dislocation|composite"),
                             "Select what type of viscosity law to use between diffusion, "
                             "dislocation and composite options. Soon there will be an option "
                             "to select a specific flow law for each assigned composition ");
          prm.declare_entry ("Yield mechanism", "drucker",
                             Patterns::Selection("drucker|limiter"),
                             "Select what type of yield mechanism to use between Drucker Prager "
                             "and stress limiter options.");

          // Diffusion creep parameters
          prm.declare_entry ("Prefactors for diffusion creep", "1.5e-15",
                             Patterns::List(Patterns::Double(0)),
                             "List of viscosity prefactors, $A$, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value. "
                             "Units: $Pa^{-n_{diffusion}} m^{n_{diffusion}/m_{diffusion}} s^{-1}$");
          prm.declare_entry ("Stress exponents for diffusion creep", "1",
                             Patterns::List(Patterns::Double(0)),
                             "List of stress exponents, $n_diffusion$, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: None");
          prm.declare_entry ("Grain size exponents for diffusion creep", "3",
                             Patterns::List(Patterns::Double(0)),
                             "List of grain size exponents, $m_diffusion$, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: None");
          prm.declare_entry ("Activation energies for diffusion creep", "375e3",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation energies, $E_a$, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: $J / mol$");
          prm.declare_entry ("Activation volumes for diffusion creep", "6e-6",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation volumes, $V_a$, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value.  Units: $m^3 / mol$");

          // Dislocation creep parameters
          prm.declare_entry ("Prefactors for dislocation creep", "1.1e-16",
                             Patterns::List(Patterns::Double(0)),
                             "List of viscosity prefactors, $A$, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value. "
                             "Units: $Pa^{-n_{dislocation}} m^{n_{dislocation}/m_{dislocation}} s^{-1}$");
          prm.declare_entry ("Stress exponents for dislocation creep", "3.5",
                             Patterns::List(Patterns::Double(0)),
                             "List of stress exponents, $n_dislocation$, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: None");
          prm.declare_entry ("Activation energies for dislocation creep", "530e3",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation energies, $E_a$, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: $J / mol$");
          prm.declare_entry ("Activation volumes for dislocation creep", "1.4e-5",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation volumes, $V_a$, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value.  Units: $m^3 / mol$");


          // Plasticity parameters
          prm.declare_entry ("Angles of internal friction", "0",
                             Patterns::List(Patterns::Double(0)),
                             "List of angles of internal friction, $\\phi$, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "For a value of zero, in 2D the von Mises criterion is retrieved. "
                             "Angles higher than 30 degrees are harder to solve numerically. Units: degrees.");
          prm.declare_entry ("Cohesions", "1e20",
                             Patterns::List(Patterns::Double(0)),
                             "List of cohesions, $C$, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "The extremely large default cohesion value (1e20 Pa) prevents the viscous stress from "
                             "exceeding the yield stress. Units: $Pa$.");

          // Stress limiter parameters
          prm.declare_entry ("Stress limiter exponents", "1.0",
                             Patterns::List(Patterns::Double(0)),
                             "List of stress limiter exponents, $n_\\text{lim}$, "
                             "for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "The default value of 1 ensures the entire stress limiter term is set to 1 "
                             "and does not affect the viscosity. Units: none.");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    ViscoPlastic<dim>::parse_parameters (ParameterHandler &prm)
    {
      //increment by one for background:
      const unsigned int n_fields = this->n_compositional_fields() + 1;

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Visco Plastic");
        {

          // Reference and minimum/maximum values
          reference_T = prm.get_double("Reference temperature");
          min_strain_rate = prm.get_double("Minimum strain rate");
          ref_strain_rate = prm.get_double("Reference strain rate");
          min_visc = prm.get_double ("Minimum viscosity");
          max_visc = prm.get_double ("Maximum viscosity");
          ref_visc = prm.get_double ("Reference viscosity");

          // Equation of state parameters
          thermal_diffusivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal diffusivities"))),
                                                                          n_fields,
                                                                          "Thermal diffusivities");
          heat_capacities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Heat capacities"))),
                                                                    n_fields,
                                                                    "Heat capacities");

          // ---- Compositional parameters
          grain_size = prm.get_double("Grain size");
          densities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Densities"))),
                                                              n_fields,
                                                              "Densities");
          thermal_expansivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal expansivities"))),
                                                                          n_fields,
                                                                          "Thermal expansivities");

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

          // Rheological parameters
          if (prm.get ("Viscous flow law") == "composite")
            viscous_flow_law = composite;
          else if (prm.get ("Viscous flow law") == "diffusion")
            viscous_flow_law = diffusion;
          else if (prm.get ("Viscous flow law") == "dislocation")
            viscous_flow_law = dislocation;
          else
            AssertThrow(false, ExcMessage("Not a valid viscous flow law"));

          // Rheological parameters
          if (prm.get ("Yield mechanism") == "drucker")
            yield_mechanism = drucker_prager;
          else if (prm.get ("Yield mechanism") == "limiter")
            yield_mechanism = stress_limiter;
          else
            AssertThrow(false, ExcMessage("Not a valid yield mechanism."));

          // Rheological parameters
          // Diffusion creep parameters (Stress exponents often but not always 1)
          prefactors_diffusion = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Prefactors for diffusion creep"))),
                                                                         n_fields,
                                                                         "Prefactors for diffusion creep");
          stress_exponents_diffusion = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Stress exponents for diffusion creep"))),
                                                                               n_fields,
                                                                               "Stress exponents for diffusion creep");
          grain_size_exponents_diffusion = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Grain size exponents for diffusion creep"))),
                                                                                   n_fields,
                                                                                   "Grain size exponents for diffusion creep");
          activation_energies_diffusion = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Activation energies for diffusion creep"))),
                                                                                  n_fields,
                                                                                  "Activation energies for diffusion creep");
          activation_volumes_diffusion = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Activation volumes for diffusion creep"))),
                                                                                 n_fields,
                                                                                 "Activation volumes for diffusion creep");
          // Dislocation creep parameters (Note the lack of grain size exponents)
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
          // Plasticity parameters
          angles_internal_friction = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Angles of internal friction"))),
                                                                             n_fields,
                                                                             "Angles of internal friction");
          cohesions = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Cohesions"))),
                                                              n_fields,
                                                              "Cohesions");
          // Stress limiter parameter
          exponents_stress_limiter  = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Stress limiter exponents"))),
                                                                              n_fields,
                                                                              "Stress limiter exponents");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::strain_rate | NonlinearDependence::compositional_fields;
      this->model_dependence.density = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(ViscoPlastic,
                                   "visco plastic",
                                   "An implementation of a visco-plastic rheology with options for "
                                   "selecting dislocation creep, diffusion creep or composite "
                                   "viscous flow laws.  Plasticity limits viscous stresses through "
                                   "a Drucker Prager yield criterion. The model is incompressible. "
                                   "Note that this material model is based heavily on the "
                                   "DiffusionDislocation (Bob Myhill) and DruckerPrager "
                                   "(Anne Glerum) material models. "
                                   "\n\n "
                                   "The viscosity for dislocation or diffusion creep is defined as "
                                   "\\[v = 0.5 * A^{-\\frac{1}{n}} * d^{-\\frac{m}{n}} * "
                                   "\\dot{\\varepsilon}_{ii}^{\\frac{1-n}{n}} * "
                                   "\\exp\\left(\\frac{E + PV}{nRT}\\right)\\] "
                                   "where $A$ is the prefactor, $n$ is the stress exponent, "
                                   "$\\dot{\\varepsilon}_{ii}$ is the square root of the deviatoric "
                                   "strain rate tensor second invariant, $d$ is grain size, "
                                   "$m$ is the grain size exponent, $E$ is activation energy, "
                                   "$V$ is activation volume, $P$ is pressure, $R$ is the gas "
                                   "exponent and $T$ is temperature. "
                                   "This form of the viscosity equation is commonly used in "
                                   "geodnynamic simulations.  See, for example, Billen and Hirth "
                                   "(2007), G3, 8, Q08012."
                                   "\n\n "
                                   "One may select to use the diffusion ($v_{diff}$; $n=1$, $m!=0$), "
                                   "dislocation ($v_{disl}$, $n>1$, $m=0$) or composite "
                                   "$\\frac{v_{diff}*v_{disl}}{v_{diff}+v_{disl}}$ equation form. "
                                   "\n\n "
                                   "Viscosity is limited through one of two different 'yielding' mechanism "
                                   "\n"
                                   "Plasticity limits viscous stress through a Drucker Prager "
                                   "yield criterion, where the yield stress in 3D is  "
                                   "$\\sigma_y = \\frac{6*C*\\cos(\\phi) + 2*P*\\sin(\\phi)} "
                                   "{\\sqrt(3)*(3+\\sin(\\phi))}$ "
                                   "and "
                                   "$\\sigma_y = C\\cos(\\phi) + P\\sin(\\phi)$ "
                                   "in 2D. Above, $C$ is cohesion and $\\phi$  is the angle of "
                                   "internal friction.  Note that the 2D form is equivalent to the "
                                   "Mohr Coloumb yield surface.  If $\\phi$ is 0, the yield stress "
                                   "is fixed and equal to the cohesion (Von Mises yield criterion) "
                                   "When the viscous stress ($2v{\\varepsilon}_{ii}$) "
                                   "the yield stress, the viscosity is rescaled back to the yield "
                                   "surface: $v_{y}=\\sigma_{y}/(2{\\varepsilon}_{ii})$. "
                                   "This form of plasticity is commonly used in geodynamic models "
                                   "See, for example, Thieulot, C. (2011), PEPI 188, pp. 47-68. "
                                   "\n\n"
                                   "Viscous stress may also be limited by a non-linear stress limiter "
                                   "that has a form similar to the Peierls creep mechanism. "
                                   "This stress limiter assigns an effective viscosity "
                                   "$\\sigma_eff = \\frac{\\tau_y}{2*\\varepsilon_y} "
                                   "{\\frac{\\varepsilon_ii}{\\varepislon_y}}^{\\frac{1}{n_y}-1}$ "
                                   "Above $\\tau_y$ is a yield stress, $\\varepsilon_y$ is the "
                                   "reference strain rate, $\\varepsilon_{ii}$ is the strain rate "
                                   "and $n_y$ is the stress limiter exponent.  The yield stress, "
                                   "$\\tau_y$, is defined through the Drucker Prager yield criterion "
                                   "formulation.  This method of limiting viscous stress has been used "
                                   "in various forms within the geodynamic literature, including "
                                   "Christensen (1992), JGR, 97(B2), pp. 2015-2036; "
                                   "Cizkova and Bina (2013), EPSL, 379, pp. 95-103; "
                                   "Cizkova and Bina (2015), EPSL, 430, pp. 408-415. "
                                   "When $n_y$ is 1, it essentially becomes a linear viscosity model,"
                                   "and in the limit $n_y\rightarrow \\infty$ it converges to the "
                                   "standard viscosity rescaling method (concretely, values $n_y>20$"
                                   "are large enough)."
                                   "\n\n "
                                   "Compositional fields can each be assigned individual values of"
                                   "thermal diffusivity, heat capacity, density, thermal "
                                   "expansivity and rheological parameters. "
                                   "\n\n "
                                   "If more than one compositional field is present at a given "
                                   "point, viscosities are averaged with an arithmetic, geometric "
                                   "harmonic (default) or maximum composition scheme. "
                                   "\n\n "
                                   "The value for the components of this formula and additional "
                                   "parameters are read from the parameter file in subsection "
                                   " 'Material model/Visco Plastic'.")
  }
}
