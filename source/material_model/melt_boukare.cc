/*
  Copyright (C) 2015 - 2019 by the authors of the ASPECT code.

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


#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/material_model/melt_boukare.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/numerics/fe_field_function.h>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    MeltBoukare<dim>::initialize()
    {
      // compute parameters for the modified Tait equation of state for the different endmembers
      // derived from the isothermal bulk modulus and its two first pressure derivatives
      // EQ 4 from Holland and Powell, 2011

      const unsigned int n_endmembers = reference_bulk_moduli.size();

      tait_parameters_a.resize(n_endmembers);
      tait_parameters_b.resize(n_endmembers);
      tait_parameters_c.resize(n_endmembers);

      for (unsigned int i=0; i<n_endmembers; ++i)
        {
          tait_parameters_a[i] = (1. + bulk_modulus_pressure_derivatives[i]) /
                                 (1. + bulk_modulus_pressure_derivatives[i] + reference_bulk_moduli[i] * bulk_modulus_second_pressure_derivatives[i]);
          tait_parameters_b[i] = bulk_modulus_pressure_derivatives[i] / reference_bulk_moduli[i]
                                 - bulk_modulus_second_pressure_derivatives[i] / (1. + bulk_modulus_pressure_derivatives[i]);
          tait_parameters_c[i] = (1. + bulk_modulus_pressure_derivatives[i] + reference_bulk_moduli[i] * bulk_modulus_second_pressure_derivatives[i]) /
                                 (std::pow(bulk_modulus_pressure_derivatives[i],2) + bulk_modulus_pressure_derivatives[i]
                                  - reference_bulk_moduli[i] * bulk_modulus_second_pressure_derivatives[i]);
        }
    }

    template <int dim>
    double
    MeltBoukare<dim>::
    reference_viscosity () const
    {
      return eta_0;
    }

    template <int dim>
    double
    MeltBoukare<dim>::
    reference_darcy_coefficient () const
    {
      // 0.01 = 1% melt
      return reference_permeability * std::pow(0.01,3.0) / eta_f;
    }

    template <int dim>
    bool
    MeltBoukare<dim>::
    is_compressible () const
    {
      return true;
    }


    template <int dim>
    double
    MeltBoukare<dim>::
    endmember_thermal_energy (const double temperature,
                              const unsigned int i) const
    {
      AssertThrow(temperature > 0.0,
                  ExcMessage("The temperature has to be larger than 0!"));

      const double energy = 3. * number_of_atoms[i] * constants::gas_constant * Einstein_temperatures[i]
                            * (0.5 + 1. / (std::exp(Einstein_temperatures[i] / temperature) - 1.0));

      return energy;
    }


    template <int dim>
    double
    MeltBoukare<dim>::
    endmember_molar_heat_capacity (const double temperature,
                                   const unsigned int i) const
    {
      AssertThrow(temperature > 0.0,
                  ExcMessage("The temperature has to be larger than 0!"));

      const double relative_T = Einstein_temperatures[i] / temperature;
      const double heat_capacity = 3. * number_of_atoms[i] * constants::gas_constant * std::pow(relative_T, 2)
                                   * std::exp(relative_T) / std::pow(std::exp(relative_T) - 1.0, 2);

      return heat_capacity;
    }


    template <int dim>
    double
    MeltBoukare<dim>::
    endmember_thermal_pressure (const double temperature,
                                const unsigned int i) const
    {
      const double thermal_energy = endmember_thermal_energy(temperature, i);
      const double heat_capacity = endmember_molar_heat_capacity(reference_temperature, i);
      const double thermal_pressure = reference_thermal_expansivities[i] * reference_bulk_moduli[i] / heat_capacity * thermal_energy;

      return thermal_pressure;
    }


    template <int dim>
    double
    MeltBoukare<dim>::
    endmember_enthalpy_thermal_addition (const double temperature,
                                         const unsigned int i) const
    {
      const double addition = reference_specific_heats[i] * temperature
                              + 0.5 * specific_heat_linear_coefficients[i] * std::pow(temperature, 2.)
                              - specific_heat_second_coefficients[i] / temperature
                              + 2. * specific_heat_third_coefficients[i] * sqrt(temperature)
                              - (reference_specific_heats[i] * reference_temperature
                                 + 0.5 * specific_heat_linear_coefficients[i] * std::pow(reference_temperature, 2.)
                                 - specific_heat_second_coefficients[i] / reference_temperature
                                 + 2.0 * specific_heat_third_coefficients[i] * sqrt(reference_temperature));

      return addition;
    }


    template <int dim>
    double
    MeltBoukare<dim>::
    endmember_entropy_thermal_addition (const double temperature,
                                        const unsigned int i) const
    {
      const double addition = reference_specific_heats[i] * std::log(temperature)
                              + specific_heat_linear_coefficients[i] * temperature
                              - 0.5 * specific_heat_second_coefficients[i] / std::pow(temperature, 2.)
                              - 2.0 * specific_heat_third_coefficients[i] / sqrt(temperature)
                              - (reference_specific_heats[i] * std::log(reference_temperature)
                                 + specific_heat_linear_coefficients[i] * reference_temperature
                                 - 0.5 * specific_heat_second_coefficients[i] / std::pow(reference_temperature, 2.)
                                 - 2.0 * specific_heat_third_coefficients[i] / sqrt(reference_temperature));

      return addition;
    }


    template <int dim>
    double
    MeltBoukare<dim>::
    melt_fraction (const double temperature,
                   const double pressure,
                   const double depletion) const
    {
      return 0.0;
    }


    template <int dim>
    void
    MeltBoukare<dim>::
    melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                    std::vector<double> &melt_fractions) const
    {
      double depletion = 0.0;

      for (unsigned int q=0; q<in.temperature.size(); ++q)
        {
          if (this->include_melt_transport())
            {
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
              const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");
              depletion = in.composition[q][peridotite_idx] - in.composition[q][porosity_idx];
            }
          melt_fractions[q] = this->melt_fraction(in.temperature[q],
                                                  std::max(0.0, in.pressure[q]),
                                                  depletion);
        }
      return;
    }


    template <int dim>
    void
    MeltBoukare<dim>::
    evaluate(const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out) const
    {
      std::vector<double> old_porosity(in.position.size());

      ReactionRateOutputs<dim> *reaction_rate_out = out.template get_additional_output<ReactionRateOutputs<dim> >();

      for (unsigned int q=0; q<in.position.size(); ++q)
        {
          if (this->include_melt_transport())
            {
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
              old_porosity[q] = in.composition[q][porosity_idx];
            }

          out.densities[q] = 0.0;
          out.thermal_expansion_coefficients[q] = 0.0;
          out.specific_heat[q] = 0.0;
          out.compressibilities[q] = 0.0;


          // TODO: more descriptive variable names and remove excess brackets once we have a unit test

          const unsigned int n_endmembers = molar_masses.size();
          for (unsigned int i=0; i<n_endmembers; ++i)
            {
              // TODO: This is just for individual endmembers, we need to average this using the composition
              const double a = tait_parameters_a[i];
              const double b = tait_parameters_b[i];
              const double c = tait_parameters_c[i];


              const double Pth = endmember_thermal_pressure(in.temperature[q], i) - endmember_thermal_pressure(reference_temperature, i);

              const double ksi_over_ksi_0 = endmember_molar_heat_capacity(in.temperature[q], i)
                                            / endmember_molar_heat_capacity(reference_temperature, i) ;

              // Integrate the gibbs free energy along the isobaric path from (P_ref, T_ref) to (P_ref, T_final)
              const double G_Pref_Tf = reference_enthalpies[i] + endmember_enthalpy_thermal_addition(in.temperature[q], i)
                                       - in.temperature[q] * (reference_entropies[i] + endmember_entropy_thermal_addition(in.temperature[q], i));

              // Integrate the gibbs free energy along the isothermal path from (P_ref, T_final) to (P_final, T_final)
              double intVdP = 0.;
              double dintVdpdT = 0.;

              if (in.pressure[q] != reference_pressure)
                {
                  intVdP = in.pressure[q] * reference_volumes[i] * (1. - a +
                                                                           (a * (std::pow((1. - b * Pth), 1. - c) -
                                                                                 std::pow((1. + b * (in.pressure[q] - Pth)), 1. - c)) /
                                                                            (b * (c - 1.) * in.pressure[q])));
                  intVdP -= reference_pressure * reference_volumes[i] * (1. - a +
                                                                         (a * (std::pow((1. - b * Pth), 1. - c) -
                                                                               std::pow((1. + b * (reference_pressure - Pth)), 1. - c)) /
                                                                          (b * (c - 1.) * reference_pressure)));

                  const double prefactor = reference_volumes[i] * reference_thermal_expansivities[i] * reference_bulk_moduli[i] * a * ksi_over_ksi_0;
                  dintVdpdT = prefactor * (std::pow(1. + b * (in.pressure[q] - Pth), -c) - std::pow(1. + b * (reference_pressure - Pth), -c));
                }


              double gibbs = G_Pref_Tf + intVdP;
              double entropy =  reference_entropies[i] + endmember_entropy_thermal_addition(in.temperature[q], i) + dintVdpdT;
              double volume = reference_volumes[i]*(1 - a * (1. - std::pow((1. + b * (in.pressure[q] - Pth)), -1.0 * c)));
              double density = molar_masses[i]/volume;

              const double isothermal_bulk_modulus = reference_bulk_moduli[i] * (1. + b * (in.pressure[q] - Pth))
                                                     * (a + (1. - a) * std::pow((1. + b * (in.pressure[q] - Pth)), c));

              const double C_V0 = endmember_molar_heat_capacity(reference_temperature, i);
              const double C_V = endmember_molar_heat_capacity(in.temperature[q], i);
              const double thermal_expansivity = reference_thermal_expansivities[i] * (C_V / C_V0) *
                                                 1. / ((1. + b * (in.pressure[q] - Pth)) *
                                                       (a + (1. - a) * std::pow((1 + b * (in.pressure[q] - Pth)), c)));

              const double Cp_ref = reference_specific_heats[i] + specific_heat_linear_coefficients[i] * in.temperature[q]
                                    + specific_heat_second_coefficients[i] * std::pow(in.temperature[q], -2.)
                                    + specific_heat_third_coefficients[i] * std::pow(in.temperature[q], -0.5);

              const double dSdT0 = reference_volumes[i] * reference_bulk_moduli[i] * std::pow(ksi_over_ksi_0 * reference_thermal_expansivities[i], 2.0)
                                   * (std::pow(1. + b * (in.pressure[q] - Pth), -1.-c) - std::pow(1. + b * (reference_pressure - Pth), -1.- c));

              const double relative_T = Einstein_temperatures[i] / in.temperature[q];
              const double dSdT = dSdT0 + dintVdpdT * (1 - 2./relative_T + 2./(std::exp(relative_T) - 1.)) * relative_T/in.temperature[q];

              const double heat_capacity_p = Cp_ref + in.temperature[q] * dSdT;



              // TODO: averaging scheme?
              // we have to get the composition index from this->introspection().compositional_index_for_name("name")
              // compositional fields do not represent endmembers, so we have to convert!
              // All properties have to be phase averaged, except for the density, which is separate for the melt and the solid
              out.densities[q] += density * in.composition[q][i];
              out.thermal_expansion_coefficients[q] += thermal_expansivity * in.composition[q][i];
              out.specific_heat[q] += heat_capacity_p * in.composition[q][i];
              out.compressibilities[q] += 1./isothermal_bulk_modulus * in.composition[q][i];
            }





          if (this->include_melt_transport() && in.strain_rate.size())
            {
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
              const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");

              // Calculate the melting rate as difference between the equilibrium melt fraction
              // and the solution of the previous time step (or the current solution, in case
              // operator splitting is used).
              // The solidus is lowered by previous melting events (fractional melting).
              const double eq_melt_fraction = melt_fraction(in.temperature[q],
                                                            this->get_adiabatic_conditions().pressure(in.position[q]),
                                                            in.composition[q][peridotite_idx] - in.composition[q][porosity_idx]);
              double porosity_change = eq_melt_fraction - old_porosity[q];

              // do not allow negative porosity
              if (old_porosity[q] + porosity_change < 0)
                porosity_change = -old_porosity[q];

              for (unsigned int c=0; c<in.composition[q].size(); ++c)
                {
                  out.reaction_terms[q][c] = 0.0;

                  // fill reaction rate outputs if the model uses operator splitting
                  if (this->get_parameters().use_operator_splitting)
                    {
                      if (reaction_rate_out != nullptr)
                        {
                          if (c == peridotite_idx && this->get_timestep_number() > 0)
                            reaction_rate_out->reaction_rates[q][c] = porosity_change / melting_time_scale
                                                                      - in.composition[q][peridotite_idx] * trace(in.strain_rate[q]);
                          else if (c == porosity_idx && this->get_timestep_number() > 0)
                            reaction_rate_out->reaction_rates[q][c] = porosity_change / melting_time_scale;
                          else
                            reaction_rate_out->reaction_rates[q][c] = 0.0;
                        }
                      out.reaction_terms[q][c] = 0.0;
                    }
                }

              const double porosity = std::min(1.0, std::max(in.composition[q][porosity_idx],0.0));
              out.viscosities[q] = eta_0 * exp(- alpha_phi * porosity);
            }
          else
            {
              out.viscosities[q] = eta_0;

              // no melting/freezing is used in the model --> set all reactions to zero
              for (unsigned int c=0; c<in.composition[q].size(); ++c)
                {
                  out.reaction_terms[q][c] = 0.0;

                  if (reaction_rate_out != nullptr)
                    reaction_rate_out->reaction_rates[q][c] = 0.0;
                }
            }

          out.entropy_derivative_pressure[q]    = 0.0;
          out.entropy_derivative_temperature[q] = 0.0;
          out.thermal_conductivities[q] = thermal_conductivity;

          double visc_temperature_dependence = 1.0;
          const double delta_temp = in.temperature[q]-this->get_adiabatic_conditions().temperature(in.position[q]);
          visc_temperature_dependence = std::max(std::min(std::exp(-thermal_viscosity_exponent*delta_temp/this->get_adiabatic_conditions().temperature(in.position[q])),1e4),1e-4);
          out.viscosities[q] *= visc_temperature_dependence;
        }

      // fill melt outputs if they exist
      MeltOutputs<dim> *melt_out = out.template get_additional_output<MeltOutputs<dim> >();

      if (melt_out != nullptr)
        {
          const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

          for (unsigned int q=0; q<in.position.size(); ++q)
            {
              double porosity = std::max(in.composition[q][porosity_idx],0.0);

              melt_out->fluid_viscosities[q] = eta_f;
              melt_out->permeabilities[q] = reference_permeability * std::pow(porosity,3) * std::pow(1.0-porosity,2);
              melt_out->fluid_density_gradients[q] = Tensor<1,dim>();

              melt_out->fluid_densities[q] = 0.0;
              melt_out->compaction_viscosities[q] = xi_0 * exp(- alpha_phi * porosity);

              const double delta_temp = in.temperature[q]-this->get_adiabatic_conditions().temperature(in.position[q]);
              const double visc_temperature_dependence = std::max(std::min(std::exp(-thermal_bulk_viscosity_exponent*delta_temp/this->get_adiabatic_conditions().temperature(in.position[q])),1e4),1e-4);

              melt_out->compaction_viscosities[q] *= visc_temperature_dependence;
            }
        }
    }



    template <int dim>
    void
    MeltBoukare<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt boukare");
        {
          prm.declare_entry ("Reference shear viscosity", "5e20",
                             Patterns::Double (0),
                             "The value of the constant viscosity $\\eta_0$ of the solid matrix. "
                             "This viscosity may be modified by both temperature and porosity "
                             "dependencies. Units: $Pa \\, s$.");
          prm.declare_entry ("Reference bulk viscosity", "1e22",
                             Patterns::Double (0),
                             "The value of the constant bulk viscosity $\\xi_0$ of the solid matrix. "
                             "This viscosity may be modified by both temperature and porosity "
                             "dependencies. Units: $Pa \\, s$.");
          prm.declare_entry ("Reference melt viscosity", "10",
                             Patterns::Double (0),
                             "The value of the constant melt viscosity $\\eta_f$. Units: $Pa \\, s$.");
          prm.declare_entry ("Exponential melt weakening factor", "27",
                             Patterns::Double (0),
                             "The porosity dependence of the viscosity. Units: dimensionless.");
          prm.declare_entry ("Thermal viscosity exponent", "0.0",
                             Patterns::Double (0),
                             "The temperature dependence of the shear viscosity. Dimensionless exponent. "
                             "See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\beta$ there.");
          prm.declare_entry ("Thermal bulk viscosity exponent", "0.0",
                             Patterns::Double (0),
                             "The temperature dependence of the bulk viscosity. Dimensionless exponent. "
                             "See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\beta$ there.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$. "
                             "Units: $W/m/K$.");
          prm.declare_entry ("Reference permeability", "1e-8",
                             Patterns::Double(),
                             "Reference permeability of the solid host rock."
                             "Units: $m^2$.");
          prm.declare_entry ("Include melting and freezing", "true",
                             Patterns::Bool (),
                             "Whether to include melting and freezing (according to a simplified "
                             "linear melting approximation in the model (if true), or not (if "
                             "false).");
          prm.declare_entry ("Melting time scale for operator splitting", "1e3",
                             Patterns::Double (0),
                             "In case the operator splitting scheme is used, the porosity field can not "
                             "be set to a new equilibrium melt fraction instantly, but the model has to "
                             "provide a melting time scale instead. This time scale defines how fast melting "
                             "happens, or more specifically, the parameter defines the time after which "
                             "the deviation of the porosity from the equilibrium melt fraction will be "
                             "reduced to a fraction of $1/e$. So if the melting time scale is small compared "
                             "to the time step size, the reaction will be so fast that the porosity is very "
                             "close to the equilibrium melt fraction after reactions are computed. Conversely, "
                             "if the melting time scale is large compared to the time step size, almost no "
                             "melting and freezing will occur."
                             "\n\n"
                             "Also note that the melting time scale has to be larger than or equal to the reaction "
                             "time step used in the operator splitting scheme, otherwise reactions can not be "
                             "computed. If the model does not use operator splitting, this parameter is not used. "
                             "Units: yr or s, depending on the ``Use years "
                             "in output instead of seconds'' parameter.");
          prm.declare_entry ("Reference temperature", "298.15",
                             Patterns::Double(),
                             "Reference temperature used to compute the material properties"
                             "of the different endmember components."
                             "Units: K.");
          prm.declare_entry ("Reference pressure", "136.e9",
                             Patterns::Double(),
                             "Reference pressure used to compute the material properties"
                             "of the different endmember components."
                             "Units: Pa.");

          prm.declare_entry ("Molar masses", "0.1003887",
                             Patterns::List(Patterns::Double(0)),
                             "Molar masses of the different endmembers"
                             "Units: kg/mol.");
          prm.declare_entry ("Number of atoms", "5.0",
                             Patterns::List(Patterns::Double(0)),
                             "Number of atoms per in the formula of each endmember."
                             "Units: none.");
          prm.declare_entry ("Reference volumes", "2.445e-05",
                             Patterns::List(Patterns::Double(0)),
                             "Reference volumes of the different endmembers."
                             "Units: $m^3$.");
          prm.declare_entry ("Reference thermal expansivities", "1.87e-05",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal expansivities for each different endmember at the reference temperature "
                             "and reference pressure."
                             "Units: 1/K.");
          prm.declare_entry ("Reference bulk moduli", "2.51e+11",
                             Patterns::List(Patterns::Double(0)),
                             "List of bulk moduli for each different endmember at the reference temperature "
                             "and reference pressure."
                             "Units: Pa.");
          prm.declare_entry ("First derivatives of the bulk modulus", "4.14",
                             Patterns::List(Patterns::Double()),
                             "The pressure derivative of the bulk modulus at the reference temperature and reference "
                             "pressure for each different endmember component."
                             "Units: none.");
          prm.declare_entry ("Second derivatives of the bulk modulus", "-1.6e-11",
                             Patterns::List(Patterns::Double()),
                             "The second pressure derivative of the bulk modulus at the reference temperature and reference "
                             "pressure for each different endmember component."
                             "Units: 1/Pa.");
          prm.declare_entry ("Einstein temperatures", "560.970464135021",
                             Patterns::List(Patterns::Double(0)),
                             "List of Einstein temperatures for each different endmember."
                             "Units: K.");
          prm.declare_entry ("Reference enthalpies", "-1442310.0",
                             Patterns::List(Patterns::Double()),
                             "List of enthalpies at the reference temperature and reference "
                             "pressure for each different endmember component."
                             "Units: J/mol.");
          prm.declare_entry ("Reference entropies", "62.6",
                             Patterns::List(Patterns::Double(0)),
                             "List of entropies at the reference temperature and reference "
                             "pressure for each different endmember component."
                             "Units: J/K/mol.");
          prm.declare_entry ("Reference specific heat capacities", "149.3",
                             Patterns::List(Patterns::Double(0)),
                             "List of specific heat capacities for each different endmember at the reference temperature "
                             "and reference pressure."
                             "Units: J/kg/K.");
          prm.declare_entry ("Linear coefficients for specific heat polynomial", "0.002918",
                             Patterns::List(Patterns::Double()),
                             "The first of three coefficients that are used to compute the specific heat capacities for each "
                             "different endmember at the reference temperature and reference pressure. "
                             "This coefficient describes the linear part of the temperature dependence. "
                             "Units: J/kg/K/K.");
          prm.declare_entry ("Second coefficients for specific heat polynomial", "-2983000.0",
                             Patterns::List(Patterns::Double()),
                             "The second of three coefficients that are used to compute the specific heat capacities for each "
                             "different endmember at the reference temperature and reference pressure. This coefficient describes "
                             "the part of the temperature dependence that scales as the inverse of the square of the temperature. "
                             "Units: J K/kg.");
          prm.declare_entry ("Third coefficients for specific heat polynomial", "-799.1",
                             Patterns::List(Patterns::Double()),
                             "The third of three coefficients that are used to compute the specific heat capacities for each "
                             "different endmember at the reference temperature and reference pressure. This coefficient describes "
                             "the part of the temperature dependence that scales as the inverse of the square root of the temperature"
                             "Units: TODO.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    MeltBoukare<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt boukare");
        {
          eta_0                             = prm.get_double ("Reference shear viscosity");
          xi_0                              = prm.get_double ("Reference bulk viscosity");
          eta_f                             = prm.get_double ("Reference melt viscosity");
          reference_permeability            = prm.get_double ("Reference permeability");
          thermal_viscosity_exponent        = prm.get_double ("Thermal viscosity exponent");
          thermal_bulk_viscosity_exponent   = prm.get_double ("Thermal bulk viscosity exponent");
          thermal_conductivity              = prm.get_double ("Thermal conductivity");
          alpha_phi                         = prm.get_double ("Exponential melt weakening factor");
          include_melting_and_freezing      = prm.get_bool ("Include melting and freezing");
          melting_time_scale                = prm.get_double ("Melting time scale for operator splitting");

          reference_temperature             = prm.get_double ("Reference temperature");
          reference_pressure                = prm.get_double ("Reference pressure");

          if (this->convert_output_to_years() == true)
            melting_time_scale *= year_in_seconds;

          AssertThrow(this->get_parameters().use_operator_splitting,
                      ExcMessage("The melt boukare material model has to be used with oprator splitting."));

          AssertThrow(melting_time_scale >= this->get_parameters().reaction_time_step,
                      ExcMessage("The reaction time step " + Utilities::to_string(this->get_parameters().reaction_time_step)
                                 + " in the operator splitting scheme is too large to compute melting rates! "
                                 "You have to choose it in such a way that it is smaller than the 'Melting time scale for "
                                 "operator splitting' chosen in the material model, which is currently "
                                 + Utilities::to_string(melting_time_scale) + "."));
          AssertThrow(melting_time_scale > 0,
                      ExcMessage("The Melting time scale for operator splitting must be larger than 0!"));

          // Equation of state parameters
          const unsigned int n_endmembers = this->include_melt_transport() ? 6 : 3;
          molar_masses = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Molar masses"))),
                                                                 n_endmembers,
                                                                 "Molar masses");
          number_of_atoms = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Number of atoms"))),
                                                                    n_endmembers,
                                                                    "Number of atoms");
          reference_volumes = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Reference volumes"))),
                                                                      n_endmembers,
                                                                      "Reference volumes");
          reference_thermal_expansivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Reference thermal expansivities"))),
                                                                                    n_endmembers,
                                                                                    "Reference thermal expansivities");
          reference_bulk_moduli = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Reference bulk moduli"))),
                                                                          n_endmembers,
                                                                          "Reference bulk moduli");
          bulk_modulus_pressure_derivatives = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("First derivatives of the bulk modulus"))),
                                                                                      n_endmembers,
                                                                                      "First derivatives of the bulk modulus");
          bulk_modulus_second_pressure_derivatives = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Second derivatives of the bulk modulus"))),
                                                     n_endmembers,
                                                     "Second derivatives of the bulk modulus");
          Einstein_temperatures = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Einstein temperatures"))),
                                                                          n_endmembers,
                                                                          "Einstein temperatures");
          reference_enthalpies = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Reference enthalpies"))),
                                                                         n_endmembers,
                                                                         "Reference enthalpies");
          reference_entropies = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Reference entropies"))),
                                                                        n_endmembers,
                                                                        "Reference entropies");
          reference_specific_heats = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Reference specific heat capacities"))),
                                                                             n_endmembers,
                                                                             "Reference specific heat capacities");
          specific_heat_linear_coefficients = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Linear coefficients for specific heat polynomial"))),
                                                                                      n_endmembers,
                                                                                      "Linear coefficients for specific heat polynomial");
          specific_heat_second_coefficients = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Second coefficients for specific heat polynomial"))),
                                                                                      n_endmembers,
                                                                                      "Second coefficients for specific heat polynomial");
          specific_heat_third_coefficients = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Third coefficients for specific heat polynomial"))),
                                                                                     n_endmembers,
                                                                                     "Third coefficients for specific heat polynomial");

          //TODO: check all lists have the correct length

          // All of these are molar fractions
          // Check that all compositional fields we need exist
          AssertThrow(this->introspection().compositional_name_exists("bridgmanite"),
                      ExcMessage("Material model melt boukare only works if there is a "
                                 "compositional field called bridgmanite."));
          AssertThrow(this->introspection().compositional_name_exists("periclase"),
                      ExcMessage("Material model melt boukare with melt transport only "
                                 "works if there is a compositional field called periclase."));
          AssertThrow(this->introspection().compositional_name_exists("iron_fraction_in_bridgmanite"),
                      ExcMessage("Material model melt boukare with melt transport only "
                                 "works if there is a compositional field called iron_fraction_in_bridgmanite."));
          AssertThrow(this->introspection().compositional_name_exists("iron_fraction_in_periclase"),
                      ExcMessage("Material model melt boukare with melt transport only "
                                 "works if there is a compositional field called iron_fraction_in_periclase."));

          if (this->include_melt_transport())
            {
              AssertThrow(this->introspection().compositional_name_exists("porosity"),
                          ExcMessage("Material model melt boukare only works if there is a "
                                     "compositional field called porosity."));
              AssertThrow(this->introspection().compositional_name_exists("magnesium_fraction_in_the_melt"),
                          ExcMessage("Material model melt boukare only works if there is a "
                                     "compositional field called magnesium_fraction_in_the_melt."));
              AssertThrow(this->introspection().compositional_name_exists("iron_fraction_in_the_melt"),
                          ExcMessage("Material model melt boukare only works if there is a "
                                     "compositional field called magnesium_fraction_in_the_melt."));
            }
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    MeltBoukare<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (this->get_parameters().use_operator_splitting
          && out.template get_additional_output<ReactionRateOutputs<dim> >() == nullptr)
        {
          const unsigned int n_points = out.viscosities.size();
          out.additional_outputs.push_back(
            std::make_shared<MaterialModel::ReactionRateOutputs<dim>> (n_points, this->n_compositional_fields()));
        }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(MeltBoukare,
                                   "melt boukare",
                                   "A material model that implements a melting model following Boukare"
                                   "et al and uses it to compute the "
                                   "material parameters required for the modelling of melt transport, "
                                   "including a source term for the porosity according to a simplified.")
  }
}
