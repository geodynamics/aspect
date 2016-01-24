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

#include <aspect/material_model/diffusion_dislocation_plasticity.h>

using namespace dealii;

namespace aspect
{
  namespace
  {
    std::vector<double>
    get_vector_double (const std::string &parameter, const unsigned int n_fields, ParameterHandler &prm)
    {
      std::vector<double> parameter_list;
      parameter_list = Utilities::string_to_double(Utilities::split_string_list(prm.get (parameter)));
      if (parameter_list.size() == 1)
        parameter_list.resize(n_fields, parameter_list[0]);

      AssertThrow(parameter_list.size() == n_fields,
                  ExcMessage("Length of "+parameter+" list must be either one, or n_compositional_fields+1"));

      return parameter_list;
    }
  }

  namespace MaterialModel
  {

    template <int dim>
    std::vector<double>
    DiffusionDislocationPlasticity<dim>::
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
    DiffusionDislocationPlasticity<dim>::
    calculate_isostrain_stress_and_stress_strain_derivative ( const std::vector<double> &volume_fractions,
                                                              const double &pressure,
                                                              const double &temperature,
                                                              const double &strain_rate_ii,
                                                              std::vector<double> &stress_ii,
                                                              std::vector<double> &stress_ii_strain_rate_ii_deriv,
                                                              std::vector<double> &stress_ii_pressure_deriv) const
    {
      // This function calculates viscosities assuming that all the compositional fields
      // experience the same strain rate (isostrain).



      // Find effective viscosities for each of the individual phases
      // Viscosities should have same number of entries as compositional fields
      std::vector<double> composition_viscosities(volume_fractions.size());
      for (unsigned int j=0; j < volume_fractions.size(); ++j)
        {
          // Power law creep equation
          // strain_rate_ii_i = A_i * stress_ii_i^{n_i} * d^{-m} \exp\left(-\frac{E_i^* + PV_i^*}{n_iRT}\right)
          // where ii indicates the square root of the second invariant and
          // i corresponds to diffusion or dislocation creep


          // For diffusion creep, viscosity is grain size dependent
          double prefactor_stress_diffusion = prefactors_diffusion[j] *
                                              std::pow(grain_size, -grain_size_exponents_diffusion[j]) *
                                              std::exp(-(activation_energies_diffusion[j] + pressure*activation_volumes_diffusion[j])/
                                                       (constants::gas_constant*temperature));

          // For dislocation creep, viscosity is grain size independent (m=0)
          double prefactor_stress_dislocation = prefactors_dislocation[j] *
                                                std::exp(-(activation_energies_dislocation[j] + pressure*activation_volumes_dislocation[j])/
                                                         (constants::gas_constant*temperature));

          double prefactor_stress_plasticity = 0.5 * std::pow(strain_rate_limit,-1/n_limit) * yield_stress[j];


          std::vector<double> L(3);
          L[0] = activation_volumes_diffusion[j]/(constants::gas_constant*temperature);
          L[1] = activation_volumes_dislocation[j]/(constants::gas_constant*temperature);
          L[2] = n_limit/(constants::gas_constant*temperature);


          // Because the ratios of the diffusion and dislocation strain rates are not known, stress is also unknown
          // We use Newton's method to find the second invariant of the stress tensor.
          // Start with the assumption that all strain is accommodated by diffusion creep:
          stress_ii[j] = strain_rate_ii/prefactor_stress_diffusion;
          double stress_ii_deriv = 1/prefactor_stress_diffusion;
          double strain_rate_residual = 2*strain_rate_residual_threshold;
          double strain_rate_deriv = 0;
          unsigned int stress_iteration = 0;
          while (std::abs(strain_rate_residual) > strain_rate_residual_threshold
                 && stress_iteration < stress_max_iteration_number)
            {
              strain_rate_residual = prefactor_stress_diffusion *
                                     std::pow(stress_ii[j], stress_exponents_diffusion[j]) +
                                     prefactor_stress_dislocation *
                                     std::pow(stress_ii[j], stress_exponents_dislocation[j]) +
                                     prefactor_stress_plasticity * std::pow(stress_ii[j],n_limit) - strain_rate_ii;

              strain_rate_deriv = stress_exponents_diffusion[j] *
                                  prefactor_stress_diffusion *
                                  std::pow(stress_ii[j], stress_exponents_diffusion[j]-1) +
                                  stress_exponents_dislocation[j] *
                                  prefactor_stress_dislocation *
                                  std::pow(stress_ii[j], stress_exponents_dislocation[j]-1) +
                                  n_limit * prefactor_stress_plasticity * std::pow(stress_ii[j],n_limit-1) ;

              stress_ii[j] -= strain_rate_residual/strain_rate_deriv;

              // initalize commom pars of the derivative to save computation time
              // to save for example on computing the expensive power function.
              std::vector<double> T_2(3);
              std::vector<double> T_1(3);
              std::vector<double> T(3);
              T_2[0] = std::pow(stress_ii[j],stress_exponents_diffusion[j]-2);
              T_1[0] = T_2[0] * stress_ii[j];
              T[0]   = T_1[0] * stress_ii[j];
              T_2[1] = std::pow(stress_ii[j],stress_exponents_dislocation[j]-2);
              T_1[1] = T_2[1] * stress_ii[j];
              T[1]   = T_1[1] * stress_ii[j];
              T_2[2] = std::pow(stress_ii[j],n_limit-2);
              T_1[2] = T_2[2] * stress_ii[j];
              T[2]   = T_1[2] * stress_ii[j];


              // compute the derivative of the stress to the strain-rate
              double sumDTn = prefactor_stress_diffusion * T[0] + prefactor_stress_dislocation * T[1] + prefactor_stress_plasticity * T[2];

              double sumDnT_1 = prefactor_stress_diffusion * stress_exponents_diffusion[j] * T_1[0] +
                                prefactor_stress_dislocation * stress_exponents_dislocation[j] * T_1[1] +
                                prefactor_stress_plasticity * n_limit * T_1[2];

              double sumDnn_1T_2 = prefactor_stress_diffusion * (stress_exponents_diffusion[j]-1) * T_2[0] +
                                   prefactor_stress_dislocation * (stress_exponents_dislocation[j]-1) * T_2[1] +
                                   prefactor_stress_plasticity * (n_limit-1) * T_2[2];

              stress_ii_strain_rate_ii_deriv[j] = stress_ii_deriv - ((1-sumDnT_1*stress_ii_strain_rate_ii_deriv[j])*sumDnT_1 - ((strain_rate_ii-sumDTn)*sumDnn_1T_2*stress_ii_strain_rate_ii_deriv[j]))/(sumDnT_1*sumDnT_1);

              // Compute the derivative of the stress to the pressure
              // We assume for now that the yieldstress is a independent of pressure.
              double sumFT       = prefactor_stress_diffusion * T[0] + prefactor_stress_dislocation * T[1] + prefactor_stress_plasticity * T[2];

              double sumFT_2n_1n = prefactor_stress_diffusion * T_2[0] * (stress_exponents_diffusion[j] - 1) * stress_exponents_diffusion[j] +
                                   prefactor_stress_dislocation * T_2[1] * (stress_exponents_dislocation[j] - 1) * stress_exponents_dislocation[j] +
                                   prefactor_stress_plasticity * T_2[2] * (n_limit -1) * n_limit;

              double sumFT_1Ln   = prefactor_stress_diffusion * T_1[0] * L[0] * stress_exponents_diffusion[j] +
                                   prefactor_stress_dislocation * T_1[1] * L[1] * stress_exponents_dislocation[j] +
                                   prefactor_stress_plasticity * T_1[2] * L[2] * n_limit;

              double sumFT_1n    = prefactor_stress_diffusion * T_1[0] * stress_exponents_diffusion[j] +
                                   prefactor_stress_dislocation * T_1[1] * stress_exponents_dislocation[j] +
                                   prefactor_stress_plasticity * T_1[2] * n_limit;

              double sumFT_2L    = prefactor_stress_diffusion * T_2[0] * L[0] +
                                   prefactor_stress_dislocation * T_2[1] * L[1] +
                                   prefactor_stress_plasticity * T_2[2] * L[2];

              stress_ii_pressure_deriv[j] = (((-sumFT + strain_rate_ii) * (sumFT_2n_1n * stress_ii_pressure_deriv[j] - sumFT_1Ln)) /
                                             (sumFT_1n*sumFT_1n))-(sumFT_1n*stress_ii_pressure_deriv[j] + sumFT_2L)/sumFT_1n;

              stress_iteration += 1;
            }


        }
      //return composition_viscosities;
    }

    template <int dim>
    void
    DiffusionDislocationPlasticity<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      //set up additional output for the derivatives
      MaterialModelDerivatives<dim> *derivatives;
      derivatives = out.template get_additional_output<MaterialModelDerivatives<dim> >();

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
              std::vector<double> stress_ii(volume_fractions.size());
              std::vector<double> stress_ii_strain_rate_ii_deriv(volume_fractions.size());
              std::vector<double> stress_ii_pressure_deriv(volume_fractions.size());

              // If strain rate is zero (like during the first time step) set it to some very small number
              // to prevent a division-by-zero, and a floating point exception.
              // Otherwise, calculate the square-root of the norm of the second invariant of the deviatoric-
              // strain rate (often simplified as epsilondot_ii)
              const double strain_rate_ii = std::max(std::sqrt(std::fabs(second_invariant(deviator(in.strain_rate[i])))),
                                                     min_strain_rate * min_strain_rate);

              // Currently, the viscosities for each of the compositional fields are calculated assuming
              // isostrain amongst all compositions, allowing calculation of the viscosity ratio.
              // TODO: This is only consistent with viscosity averaging if the arithmetic averaging
              // scheme is chosen. It would be useful to have a function to calculate isostress viscosities.
              calculate_isostrain_stress_and_stress_strain_derivative(volume_fractions, pressure, temperature, strain_rate_ii, stress_ii, stress_ii_strain_rate_ii_deriv,stress_ii_pressure_deriv);

              std::vector<double> viscosity_strain_rate_ii_deriv(volume_fractions.size());
              std::vector<double> viscosity_pressure_deriv(volume_fractions.size());
              // Compute the viscosity
              std::vector<double> viscosities(volume_fractions.size());
              for (unsigned int j=0; j < volume_fractions.size(); ++j)
                {
                  // The effective viscosity, with minimum and maximum bounds
                  viscosities[j] = volume_fractions[j] * std::pow(std::min(std::max(stress_ii[j]/strain_rate_ii/2, min_visc), max_visc),composition_averaging_exponent);
                  out.viscosities[i] += viscosities[j];

                  // For diffusion creep, viscosity is grain size dependent
                  double prefactor_stress_diffusion = prefactors_diffusion[j] *
                                                      std::pow(grain_size, -grain_size_exponents_diffusion[j]) *
                                                      std::exp(-(activation_energies_diffusion[j] + pressure*activation_volumes_diffusion[j])/
                                                               (constants::gas_constant*temperature));

                  // For dislocation creep, viscosity is grain size independent (m=0)
                  double prefactor_stress_dislocation = prefactors_dislocation[j] *
                                                        std::exp(-(activation_energies_dislocation[j] + pressure*activation_volumes_dislocation[j])/
                                                                 (constants::gas_constant*temperature));

                  double prefactor_stress_plasticity = 0.5 * std::pow(strain_rate_limit,-1/n_limit) * yield_stress[j];

                  double sumFT1_n = prefactor_stress_diffusion * std::pow(stress_ii[j],stress_exponents_diffusion[j]) * (1-stress_exponents_diffusion[j]) +
                                    prefactor_stress_dislocation * std::pow(stress_ii[j],stress_exponents_dislocation[j]) * (1-stress_exponents_dislocation[j]) +
                                    prefactor_stress_plasticity * std::pow(stress_ii[j],n_limit) * (1-n_limit);

                  double sumFT = prefactor_stress_diffusion * std::pow(stress_ii[j],stress_exponents_diffusion[j]) +
                                 prefactor_stress_dislocation * std::pow(stress_ii[j],stress_exponents_dislocation[j]) +
                                 prefactor_stress_plasticity * std::pow(stress_ii[j],n_limit);

                  double sumFn_1 = prefactor_stress_diffusion * (stress_exponents_diffusion[j]-1) +
                                   prefactor_stress_dislocation * (stress_exponents_dislocation[j]-1) +
                                   prefactor_stress_plasticity *  (n_limit-1);

                  double sumFT1 = prefactor_stress_diffusion * std::pow(stress_ii[j],stress_exponents_diffusion[j]+1) +
                                  prefactor_stress_dislocation * std::pow(stress_ii[j],stress_exponents_dislocation[j]+1) +
                                  prefactor_stress_plasticity * std::pow(stress_ii[j],n_limit+1);

                  viscosity_strain_rate_ii_deriv[j] = (sumFT1_n/(sumFT*sumFT))*stress_ii_strain_rate_ii_deriv[j];

                  viscosity_pressure_deriv[j] = (sumFn_1 * constants::gas_constant * temperature * stress_ii_pressure_deriv[j] + sumFT1) /
                                                (2 * constants::gas_constant * temperature * sumFT * sumFT);


                  // TODO: review if this this old description is still valid when plasticity is added.
                  // The isostrain condition implies that the viscosity averaging should be arithmetic (see above).
                  // We have given the user freedom to apply alternative bounds, because in diffusion-dominated
                  // creep (where n_diff=1) viscosities are stress and strain-rate independent, so the calculation
                  // of compositional field viscosities is consistent with any averaging scheme.
                  //out.viscosities[i] = average_value(composition, composition_viscosities, viscosity_averaging);
                }


              out.viscosities[i] = 1/volume_fractions.size() * pow(out.viscosities[i],composition_averaging_exponent_inverse);

              // Compute the viscosity derivative
              derivatives->dviscosities_dvelocity[i] = 0.0;
              double sum_c_eta_p = 0;
              double sum_c_p_detadepsilon_p_1 = 0;
              double sum_c_p_detadP_p_1 = 0;
              for (unsigned int j=0; j < volume_fractions.size(); ++j)
                {
                  //TODO: look at this part again!
                  sum_c_eta_p     += volume_fractions[j] * std::pow(viscosities[j], composition_averaging_exponent);

                  sum_c_p_detadepsilon_p_1 += volume_fractions[j] * composition_averaging_exponent * std::pow(viscosity_strain_rate_ii_deriv[j], composition_averaging_exponent-1);
                  sum_c_p_detadP_p_1 += volume_fractions[j] * composition_averaging_exponent * std::pow(viscosity_pressure_deriv[j], composition_averaging_exponent-1);
                }
              derivatives->dviscosities_dstrain_rate[i] = 1/volume_fractions.size() * composition_averaging_exponent_inverse *
                                                          pow(sum_c_eta_p,1/composition_averaging_exponent) * sum_c_p_detadepsilon_p_1;

              derivatives->dviscosities_dpressure[i] = 1/volume_fractions.size() * composition_averaging_exponent_inverse *
                                                       pow(sum_c_eta_p,1/composition_averaging_exponent) * sum_c_p_detadP_p_1;

              derivatives->dviscosities_dtemperature[i] = 0.0;
              //derivatives->dviscosities_dcompositions[i] = std::vector<double>(in.composition[i].size(),0);
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
    DiffusionDislocationPlasticity<dim>::
    reference_viscosity () const
    {
      return ref_visc;
    }

    template <int dim>
    double
    DiffusionDislocationPlasticity<dim>::
    reference_density () const
    {
      return densities[0];
    }

    template <int dim>
    bool
    DiffusionDislocationPlasticity<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    DiffusionDislocationPlasticity<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Diffusion dislocation plasticity");
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

          // Viscosity iteration parameters
          prm.declare_entry ("Strain rate residual tolerance", "1e-22", Patterns::Double(0),
                             "Tolerance for correct diffusion/dislocation strain rate ratio.");
          prm.declare_entry ("Maximum strain rate ratio iterations", "40", Patterns::Integer(0),
                             "Maximum number of iterations to find the correct "
                             "diffusion/dislocation strain rate ratio.");

          // Equation of state parameters
          prm.declare_entry ("Thermal diffusivity", "0.8e-6", Patterns::Double(0), "Units: $m^2/s$");
          prm.declare_entry ("Heat capacity", "1.25e3", Patterns::Double(0), "Units: $J / (K * kg)$");
          prm.declare_entry ("Densities", "3300.",
                             Patterns::List(Patterns::Double(0)),
                             "List of densities, $\\rho$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: $kg / m^3$");
          prm.declare_entry ("Thermal expansivities", "3.5e-5",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal expansivities for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: $1 / K$");

          // Rheological parameters
          prm.declare_entry ("Grain size", "1e-3", Patterns::Double(0), "Units: $m$");
          prm.declare_entry ("Composition averaging exponent", "-1",
                             Patterns::Double(),
                             "When more than one compositional field is present at a point "
                             "with different viscosities, we need to come up with an average "
                             "viscosity at that point.  Select a weighted harmonic (-1), arithmetic (1), "
                             "geometric, or maximum (very large) composition by setting the "
                             "appropeate value of the implemented general mean.");

          // plasticity parameters
          prm.declare_entry ("n limit", "50", Patterns::Integer(0),
                             "Power-law index defining the 'brittleness' of the material.");
          prm.declare_entry ("Strain rate limit", "1e-15", Patterns::Double(0),
                             "Strain rate limit");
          prm.declare_entry ("Yield stress","1e9",Patterns::List(Patterns::Double(0)),
                             "List of Yield stresses per material.");
          // Diffusion creep parameters
          prm.declare_entry ("Prefactors for diffusion creep", "1.5e-15",
                             Patterns::List(Patterns::Double(0)),
                             "List of viscosity prefactors, $A$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value. "
                             "Units: $Pa^{-n_{diffusion}} m^{n_{diffusion}/m_{diffusion}} s^{-1}$");
          prm.declare_entry ("Stress exponents for diffusion creep", "1",
                             Patterns::List(Patterns::Double(0)),
                             "List of stress exponents, $n_diffusion$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: None");
          prm.declare_entry ("Grain size exponents for diffusion creep", "3",
                             Patterns::List(Patterns::Double(0)),
                             "List of grain size exponents, $m_diffusion$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: None");
          prm.declare_entry ("Activation energies for diffusion creep", "375e3",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation energies, $E_a$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: $J / mol$");
          prm.declare_entry ("Activation volumes for diffusion creep", "6e-6",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation volumes, $V_a$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value.  Units: $m^3 / mol$");

          // Dislocation creep parameters
          prm.declare_entry ("Prefactors for dislocation creep", "1.1e-16",
                             Patterns::List(Patterns::Double(0)),
                             "List of viscosity prefactors, $A$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value. "
                             "Units: $Pa^{-n_{dislocation}} m^{n_{dislocation}/m_{dislocation}} s^{-1}$");
          prm.declare_entry ("Stress exponents for dislocation creep", "3.5",
                             Patterns::List(Patterns::Double(0)),
                             "List of stress exponents, $n_dislocation$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: None");
          prm.declare_entry ("Activation energies for dislocation creep", "530e3",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation energies, $E_a$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one values is given, then all use the same value.  Units: $J / mol$");
          prm.declare_entry ("Activation volumes for dislocation creep", "1.4e-5",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation volumes, $V_a$, for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value.  Units: $m^3 / mol$");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    DiffusionDislocationPlasticity<dim>::parse_parameters (ParameterHandler &prm)
    {
      //increment by one for background:
      const unsigned int n_fields = this->n_compositional_fields() + 1;

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Diffusion dislocation plasticity");
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

          // Iteration parameters
          strain_rate_residual_threshold = prm.get_double ("Strain rate residual tolerance");
          stress_max_iteration_number = prm.get_integer ("Maximum strain rate ratio iterations");

          // Equation of state parameters
          thermal_diffusivity = prm.get_double("Thermal diffusivity");
          heat_capacity = prm.get_double("Heat capacity");

          // ---- Compositional parameters
          grain_size = prm.get_double("Grain size");
          densities = get_vector_double("Densities", n_fields, prm);
          thermal_expansivities = get_vector_double("Thermal expansivities", n_fields, prm);

          // Rheological parameters
          composition_averaging_exponent = prm.get_double ("composition_averaging_exponent");

          // Rheological parameters
          // Plasticity parameters
          n_limit = prm.get_integer("n limit");
          strain_rate_limit = prm.get_double("Strain rate limit");
          yield_stress = get_vector_double("Yield stress", n_fields, prm);
          // Diffusion creep parameters (Stress exponents often but not always 1)
          prefactors_diffusion = get_vector_double("Prefactors for diffusion creep", n_fields, prm);
          stress_exponents_diffusion = get_vector_double("Stress exponents for diffusion creep", n_fields, prm);
          grain_size_exponents_diffusion = get_vector_double("Grain size exponents for diffusion creep", n_fields, prm);
          activation_energies_diffusion = get_vector_double("Activation energies for diffusion creep", n_fields, prm);
          activation_volumes_diffusion = get_vector_double("Activation volumes for diffusion creep", n_fields, prm);
          // Dislocation creep parameters (Note the lack of grain size exponents)
          prefactors_dislocation = get_vector_double("Prefactors for dislocation creep", n_fields, prm);
          stress_exponents_dislocation = get_vector_double("Stress exponents for dislocation creep", n_fields, prm);
          activation_energies_dislocation = get_vector_double("Activation energies for dislocation creep", n_fields, prm);
          activation_volumes_dislocation = get_vector_double("Activation volumes for dislocation creep", n_fields, prm);

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
    ASPECT_REGISTER_MATERIAL_MODEL(DiffusionDislocationPlasticity,
                                   "diffusion dislocation plasticity",
                                   " An implementation of a viscous rheology including diffusion"
                                   " and dislocation creep."
                                   " Compositional fields can each be assigned individual"
                                   " activation energies, reference densities, thermal expansivities,"
                                   " and stress exponents. The effective viscosity is defined as"
                                   " \n\n"
                                   " \\[v_\\text{eff} = \\left(\\frac{1}{v_\\text{eff}^\\text{diff}}+"
                                   " \\frac{1}{v_\\text{eff}^\\text{dis}}\\right)^{-1}\\]"
                                   " where"
                                   " \\[v_\\text{i} = 0.5 * A^{-\\frac{1}{n_i}} d^\\frac{m_i}{n_i}"
                                   " \\dot{\\varepsilon_i}^{\\frac{1-n_i}{n_i}}"
                                   " \\exp\\left(\\frac{E_i^* + PV_i^*}{n_iRT}\\right)\\]"
                                   " \n\n"
                                   " where $d$ is grain size, $i$ corresponds to diffusion or dislocation creep,"
                                   " $\\dot{\\varepsilon}$ is the square root of the second invariant of the"
                                   " strain rate tensor, $R$ is the gas constant, $T$ is temperature, "
                                   " and $P$ is pressure."
                                   " $A_i$ are prefactors, $n_i$ and $m_i$ are stress and grain size exponents"
                                   " $E_i$ are the activation energies and $V_i$ are the activation volumes."
                                   " \n\n"
                                   " The ratio of diffusion to dislocation strain rate is found by Newton's"
                                   " method, iterating to find the stress which satisfies the above equations."
                                   " The value for the components of this formula and additional"
                                   " parameters are read from the parameter file in subsection"
                                   " 'Material model/DiffusionDislocation'.")
  }
}
