/*
  Copyright (C) 2015 - 2024 by the authors of the ASPECT code.

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
#include <aspect/gravity_model/interface.h>
#include <aspect/simulator.h>

#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace
    {
      std::vector<std::string> make_boukare_additional_outputs_names()
      {
        std::vector<std::string> names;
        names.emplace_back("bulk_composition");
        names.emplace_back("molar_volatiles_in_melt");
        return names;
      }
    }



    template <int dim>
    BoukareOutputs<dim>::BoukareOutputs (const unsigned int n_points)
      :
      NamedAdditionalMaterialOutputs<dim>(make_boukare_additional_outputs_names()),
      bulk_composition(n_points, numbers::signaling_nan<double>()),
      molar_volatiles_in_melt(n_points, numbers::signaling_nan<double>())
    {}



    template <int dim>
    std::vector<double>
    BoukareOutputs<dim>::get_nth_output(const unsigned int idx) const
    {
      AssertIndexRange (idx, 2);
      if (idx == 0)
        return bulk_composition;
      else
        return molar_volatiles_in_melt;
    }



    template <int dim>
    void
    MeltBoukare<dim>::initialize()
    {
      CitationInfo::add("boukaremelt");

      // Compute parameters for the modified Tait equation of state for the different endmembers
      // derived from the isothermal bulk modulus and its two first pressure derivatives.
      // This corresponds to Equation 4 from Holland and Powell, 2011 (https://doi.org/10.1111/j.1525-1314.2010.00923.x).

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
                                 (Utilities::fixed_power<2>(bulk_modulus_pressure_derivatives[i]) + bulk_modulus_pressure_derivatives[i]
                                  - reference_bulk_moduli[i] * bulk_modulus_second_pressure_derivatives[i]);
        }

      febdg_idx = std::distance(endmember_names.begin(), find(endmember_names.begin(), endmember_names.end(), "FeSiO3_bridgmanite"));
      mgbdg_idx = std::distance(endmember_names.begin(), find(endmember_names.begin(), endmember_names.end(), "MgSiO3_bridgmanite"));
      wus_idx = std::distance(endmember_names.begin(), find(endmember_names.begin(), endmember_names.end(), "FeO_periclase"));
      per_idx = std::distance(endmember_names.begin(), find(endmember_names.begin(), endmember_names.end(), "MgO_periclase"));
      femelt_idx = std::distance(endmember_names.begin(), find(endmember_names.begin(), endmember_names.end(), "FeO_melt"));
      mgmelt_idx = std::distance(endmember_names.begin(), find(endmember_names.begin(), endmember_names.end(), "MgO_melt"));
    }



    template <int dim>
    double
    MeltBoukare<dim>::
    reference_darcy_coefficient () const
    {
      // 0.01 = 1% melt
      return reference_permeability * Utilities::fixed_power<3>(0.01) / eta_f;
    }

    template <int dim>
    bool
    MeltBoukare<dim>::
    is_compressible () const
    {
      return true;
    }


    template <int dim>
    MeltBoukare<dim>::
    EndmemberProperties::EndmemberProperties(const unsigned int n_endmembers)
      :
      volumes(n_endmembers, numbers::signaling_nan<double>()),
      gibbs_energies(n_endmembers, numbers::signaling_nan<double>()),
      entropies(n_endmembers, numbers::signaling_nan<double>()),
      thermal_expansivities(n_endmembers, numbers::signaling_nan<double>()),
      bulk_moduli(n_endmembers, numbers::signaling_nan<double>()),
      heat_capacities(n_endmembers, numbers::signaling_nan<double>())
    {}



    template <int dim>
    void
    MeltBoukare<dim>::
    fill_endmember_properties (const typename Interface<dim>::MaterialModelInputs &in,
                               const unsigned int q,
                               EndmemberProperties &properties) const
    {
      const double n_endmembers = properties.volumes.size();

      for (unsigned int i=0; i<n_endmembers; ++i)
        {
          const double a = tait_parameters_a[i];
          const double b = tait_parameters_b[i];
          const double c = tait_parameters_c[i];

          const double Pth = endmember_thermal_pressure(in.temperature[q], i) - endmember_thermal_pressure(reference_temperature, i);

          const double heat_capacity_ratio = endmember_molar_heat_capacity(in.temperature[q], i)
                                             / endmember_molar_heat_capacity(reference_temperature, i) ;

          // Integrate the gibbs free energy along the isobaric path from (P_ref, T_ref) to (P_ref, T_final)
          const double G_Pref_Tf = reference_enthalpies[i] + endmember_enthalpy_thermal_addition(in.temperature[q], i)
                                   - in.temperature[q] * (reference_entropies[i] + endmember_entropy_thermal_addition(in.temperature[q], i));

          // Integrate the gibbs free energy along the isothermal path from (P_ref, T_final) to (P_final, T_final)
          long double intVdP = 0.;
          long double dintVdpdT = 0.;

          const double pressure = this->get_adiabatic_conditions().pressure(in.position[q]);

          if (pressure != reference_pressure && pressure > 0.0)
            {
              intVdP = reference_volumes[i]
                       * ((pressure - reference_pressure) * (1. - a)
                          + (a * (std::pow((1. + b * (reference_pressure - Pth)), 1. - c) - std::pow((1. + b * (pressure - Pth)), 1. - c)) / (b * (c - 1.))));

              const double prefactor = reference_volumes[i] * reference_thermal_expansivities[i] * reference_bulk_moduli[i] * a * heat_capacity_ratio;
              dintVdpdT = prefactor * (std::pow(1. + b * (pressure - Pth), -c) - std::pow(1. + b * (reference_pressure - Pth), -c));
            }

          properties.gibbs_energies[i] = G_Pref_Tf + intVdP;
          properties.entropies[i] =  reference_entropies[i] + endmember_entropy_thermal_addition(in.temperature[q], i) + dintVdpdT;
          properties.volumes[i] = reference_volumes[i] * (1 - a * (1. - std::pow(1. + b * (pressure - Pth), -c)));

          properties.bulk_moduli[i] = reference_bulk_moduli[i] * (1. + b * (pressure - Pth))
                                      * (a + (1. - a) * std::pow(1. + b * (pressure - Pth), c));

          const double C_V0 = endmember_molar_heat_capacity(reference_temperature, i);
          const double C_V = endmember_molar_heat_capacity(in.temperature[q], i);
          properties.thermal_expansivities[i] = reference_thermal_expansivities[i] * (C_V / C_V0) *
                                                1. / ((1. + b * (pressure - Pth)) *
                                                      (a + (1. - a) * std::pow(1 + b * (pressure - Pth), c)));

          const double Cp_ref = reference_specific_heats[i] + specific_heat_linear_coefficients[i] * in.temperature[q]
                                + specific_heat_second_coefficients[i] * Utilities::fixed_power<-2>(in.temperature[q])
                                + specific_heat_third_coefficients[i] * std::pow(in.temperature[q], -0.5);

          const long double dSdT0 = reference_volumes[i] * reference_bulk_moduli[i] * Utilities::fixed_power<2>(heat_capacity_ratio * reference_thermal_expansivities[i])
                                    * (std::pow(1. + b * (pressure - Pth), -1.-c) - std::pow(1. + b * (reference_pressure - Pth), -1.- c));

          const double relative_T = Einstein_temperatures[i] / in.temperature[q];
          const double dSdT = dSdT0 + dintVdpdT * (1 - 2./relative_T + 2./(std::exp(relative_T) - 1.)) * relative_T/in.temperature[q];

          properties.heat_capacities[i] = Cp_ref + in.temperature[q] * dSdT;
        }
    }



    template <int dim>
    double
    MeltBoukare<dim>::
    endmember_thermal_energy (const double temperature,
                              const unsigned int endmember_index) const
    {
      AssertThrow(temperature > 0.0,
                  ExcMessage("The temperature has to be larger than 0, but it is "
                             + std::to_string(temperature) + " for endmember " + std::to_string(endmember_index) + "."));

      const double relative_T = Einstein_temperatures[endmember_index] / temperature;
      const double energy = 3. * number_of_atoms[endmember_index] * constants::gas_constant * Einstein_temperatures[endmember_index]
                            * (0.5 + 1. / (std::exp(relative_T) - 1.0));

      return energy;
    }



    template <int dim>
    double
    MeltBoukare<dim>::
    endmember_molar_heat_capacity (const double temperature,
                                   const unsigned int endmember_index) const
    {
      AssertThrow(temperature > 0.0,
                  ExcMessage("The temperature has to be larger than 0!"));

      const double relative_T = Einstein_temperatures[endmember_index] / temperature;
      const double heat_capacity = 3. * number_of_atoms[endmember_index] * constants::gas_constant * Utilities::fixed_power<2>(relative_T)
                                   * std::exp(relative_T) / Utilities::fixed_power<2>(std::exp(relative_T) - 1.0);

      return heat_capacity;
    }



    template <int dim>
    double
    MeltBoukare<dim>::
    endmember_thermal_pressure (const double temperature,
                                const unsigned int endmember_index) const
    {
      const double thermal_energy = endmember_thermal_energy(temperature, endmember_index);
      const double heat_capacity = endmember_molar_heat_capacity(reference_temperature, endmember_index);
      const double thermal_pressure = reference_thermal_expansivities[endmember_index] * reference_bulk_moduli[endmember_index] / heat_capacity * thermal_energy;

      return thermal_pressure;
    }



    template <int dim>
    double
    MeltBoukare<dim>::
    endmember_enthalpy_thermal_addition (const double temperature,
                                         const unsigned int i) const
    {
      const double addition = reference_specific_heats[i] * temperature
                              + 0.5 * specific_heat_linear_coefficients[i] * Utilities::fixed_power<2>(temperature)
                              - specific_heat_second_coefficients[i] / temperature
                              + 2. * specific_heat_third_coefficients[i] * std::sqrt(temperature)
                              - (reference_specific_heats[i] * reference_temperature
                                 + 0.5 * specific_heat_linear_coefficients[i] * Utilities::fixed_power<2>(reference_temperature)
                                 - specific_heat_second_coefficients[i] / reference_temperature
                                 + 2.0 * specific_heat_third_coefficients[i] * std::sqrt(reference_temperature));

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
                              - 0.5 * specific_heat_second_coefficients[i] / Utilities::fixed_power<2>(temperature)
                              - 2.0 * specific_heat_third_coefficients[i] / std::sqrt(temperature)
                              - (reference_specific_heats[i] * std::log(reference_temperature)
                                 + specific_heat_linear_coefficients[i] * reference_temperature
                                 - 0.5 * specific_heat_second_coefficients[i] / Utilities::fixed_power<2>(reference_temperature)
                                 - 2.0 * specific_heat_third_coefficients[i] / std::sqrt(reference_temperature));

      return addition;
    }



    template <int dim>
    void
    MeltBoukare<dim>::
    convert_composition_to_fraction_of_endmembers (const double temperature,
                                                   const double molar_Fe_in_solid,
                                                   const double molar_Fe_in_melt,
                                                   const std::vector<double> &endmember_gibbs_energies,
                                                   std::vector<double> &endmember_mole_fractions_per_phase,
                                                   double &molar_bridgmanite_in_solid) const
    {
      const double x_FeO = molar_Fe_in_solid * molar_FeO_in_Fe_mantle_endmember;
      const double x_MgO = (1. - molar_Fe_in_solid) * molar_MgO_in_Mg_mantle_endmember;
      const double x_SiO2 = molar_Fe_in_solid * molar_SiO2_in_Fe_mantle_endmember + (1. - molar_Fe_in_solid) * molar_SiO2_in_Mg_mantle_endmember;

      const double molar_fraction_FeO = x_FeO/(x_MgO + x_FeO);
      const double molar_fraction_SiO2 = x_SiO2/(x_MgO + x_FeO);
      molar_bridgmanite_in_solid = molar_fraction_SiO2;

      const double gibbs_energy_of_reaction = endmember_gibbs_energies[mgbdg_idx] + endmember_gibbs_energies[wus_idx]
                                              - endmember_gibbs_energies[febdg_idx] - endmember_gibbs_energies[per_idx];
      const double partition_coefficient = std::exp(gibbs_energy_of_reaction / (constants::gas_constant * temperature));

      // Solving equation 6 in Nakajima et al., 2012 for X_Fe_fp and X_Fe_pv
      // Solved using the definition of the distribution coefficient to define X_Fe_fp as a function of X_Fe_pv

      const double num_to_sqrt = -4. * molar_fraction_FeO * (partition_coefficient - 1.) * partition_coefficient * molar_fraction_SiO2
                                 + Utilities::fixed_power<2>(1. + (molar_fraction_FeO + molar_fraction_SiO2) * (partition_coefficient - 1.0));

      endmember_mole_fractions_per_phase[febdg_idx] = (-1. + molar_fraction_FeO - (molar_fraction_FeO * partition_coefficient) + molar_fraction_SiO2 - (molar_fraction_SiO2 * partition_coefficient) + std::sqrt(num_to_sqrt)) /
                                                      (2. * molar_fraction_SiO2 * (1. - partition_coefficient));

      endmember_mole_fractions_per_phase[wus_idx] = endmember_mole_fractions_per_phase[febdg_idx] / (((1. - endmember_mole_fractions_per_phase[febdg_idx]) * partition_coefficient) + endmember_mole_fractions_per_phase[febdg_idx]);

      endmember_mole_fractions_per_phase[mgbdg_idx] = 1.0 - endmember_mole_fractions_per_phase[febdg_idx];
      endmember_mole_fractions_per_phase[per_idx] = 1.0 - endmember_mole_fractions_per_phase[wus_idx];

      endmember_mole_fractions_per_phase[mgmelt_idx] = 1. - molar_Fe_in_melt;
      endmember_mole_fractions_per_phase[femelt_idx] = molar_Fe_in_melt;

      return;
    }



    template <int dim>
    double
    MeltBoukare<dim>::
    compute_melt_molar_fraction (const double porosity,
                                 const double bridgmanite_molar_fraction_in_solid,
                                 EndmemberProperties &endmembers,
                                 const std::vector<double> &endmember_mole_fractions_per_phase) const
    {
      double melt_molar_volume = 0.0;
      double solid_molar_volume = 0.0;
      for (unsigned int i=0; i<endmember_names.size(); ++i)
        {
          // We assume here that the first 2 endmembers are the (solid) bridgmanite endmembers,
          // and the 3rd and 4th are the (solid) ferropericlase endmembers. The endmember_phase_fraction_in_solid
          // variable therefore indicates the fraction of bridgmanite or ferropericlase in the solid
          // (depending on which phase the endmember belongs to).
          const double endmember_phase_fraction_in_solid = i<2 ? bridgmanite_molar_fraction_in_solid : 1.0 - bridgmanite_molar_fraction_in_solid;

          if (endmember_states[i] == EndmemberState::melt)
            melt_molar_volume += endmember_mole_fractions_per_phase[i] * endmembers.volumes[i];
          else if (endmember_states[i] == EndmemberState::solid)
            solid_molar_volume += endmember_phase_fraction_in_solid * endmember_mole_fractions_per_phase[i] * endmembers.volumes[i];
          else
            AssertThrow (false, ExcNotImplemented());
        }

      // The porosity needs to be between 0 and 1,
      // but that may not always be the case, because it is advected using the finite element.
      const double bounded_porosity = std::min(1.0, std::max(0.0, porosity));

      double melt_molar_fraction = 0.0;
      if (melt_molar_volume > 0.0)
        {
          const double n_moles_total = bounded_porosity / melt_molar_volume + (1.0 - bounded_porosity) / solid_molar_volume;
          melt_molar_fraction = bounded_porosity / (melt_molar_volume * n_moles_total);
        }

      return melt_molar_fraction;
    }



    template <int dim>
    double
    MeltBoukare<dim>::
    assert_update_is_within_0_and_1 (const double old_value,
                                     const double change_of_value) const
    {
      if (old_value + change_of_value < 0.0)
        {
          AssertThrow (false, ExcMessage("Update below 0. Proposed update: " + std::to_string(old_value + change_of_value)));
          return -old_value;
        }
      else if (old_value + change_of_value > 1.0)
        {
          AssertThrow (false, ExcMessage("Update above 1. Proposed update: " + std::to_string(old_value + change_of_value)));
          return 1.0 - old_value;
        }
      else
        return change_of_value;
    }



    template <int dim>
    double
    MeltBoukare<dim>::
    melt_fraction (const double temperature,
                   const double pressure,
                   const double molar_composition_of_bulk,
                   double &molar_volatiles_in_melt,
                   double &new_molar_composition_of_solid,
                   double &new_molar_composition_of_melt) const
    {
      if (temperature == 0.0)
        return 0.0;

      const double molar_volatiles_in_bulk = 1.e-4;
      {
        // after Phipps Morgan, Jason. "Thermodynamics of pressure release melting of a veined plum pudding mantle."
        // Geochemistry, Geophysics, Geosystems 2.4 (2001).
        // See also Appendix B of Dannberg et al., 2021. "The morphology, evolution and seismic visibility of partial
        // melt at the core-mantle boundary: Implications for ULVZs".
        const double P = pressure;                   // pressure in Pa
        const double T = temperature;                // temperature in K
        const double R = constants::gas_constant;    // Ideal Gas Constant

        // Free Energy Change Delta_G due to melting as a function of temperature and pressure
        const double dG_Fe_mantle = (Fe_mantle_melting_temperature - T) * Fe_mantle_melting_entropy
                                    + (P - melting_reference_pressure) * Fe_mantle_melting_volume;
        const double dG_Mg_mantle = (Mg_mantle_melting_temperature - T) * Mg_mantle_melting_entropy
                                    + (P - melting_reference_pressure) * Mg_mantle_melting_volume;

        // Equations (B.13) and (B.14) in Appendix B of Dannberg et al., 2021.
        const double c_Fe_endmember = std::exp(dG_Fe_mantle/(Fe_number_of_moles*R*T));
        const double c_Mg_endmember = std::exp(dG_Mg_mantle/(Mg_number_of_moles*R*T));

        // Mole composition of the solid and liquid (corresponds to the molar fraction of X_Fe, the iron-bearing endmember).
        // In addition to the Phipps Morgan model, we also include volatiles here.
        // Equations (B.8) and (B.10) to (B.12) in Appendix B of Dannberg et al., 2021.
        const double q0 = (c_Fe_endmember - c_Mg_endmember) * c_Fe_endmember/c_Mg_endmember;
        const double q1 = c_Fe_endmember * (1. - 1./c_Mg_endmember) - molar_composition_of_bulk * (c_Fe_endmember - c_Mg_endmember)/c_Mg_endmember + molar_volatiles_in_bulk * (1. - c_Fe_endmember);
        const double q2 = -molar_composition_of_bulk * (1. - 1./c_Mg_endmember);
        const double molar_composition_of_melt = (-q1 + std::sqrt(q1*q1 - 4.*q0*q2))/(2.*q0);

        const double T_Fe_mantle = dG_Fe_mantle / Fe_mantle_melting_entropy;
        const double T_Mg_mantle = dG_Mg_mantle / Mg_mantle_melting_entropy;

        double melt_molar_fraction;

        if (molar_composition_of_bulk < std::numeric_limits<double>::min())
          {
            melt_molar_fraction = T_Mg_mantle < 0.0 ? 1.0 : molar_volatiles_in_bulk;
            new_molar_composition_of_melt = 0.0;
            new_molar_composition_of_solid = 0.0;
          }
        else if (molar_composition_of_melt <= molar_composition_of_bulk
                 && std::min(T_Fe_mantle, T_Mg_mantle) < 0.0) // above the liquidus
          {
            melt_molar_fraction = 1.0;
            new_molar_composition_of_melt = molar_composition_of_bulk;
            new_molar_composition_of_solid = molar_composition_of_bulk;

            molar_volatiles_in_melt = molar_volatiles_in_bulk;
          }
        else
          {
            // Equations (B.9) in Appendix B of Dannberg et al., 2021.
            const double molar_composition_of_solid = molar_composition_of_melt * c_Fe_endmember;

            new_molar_composition_of_solid = std::max(std::min(molar_composition_of_solid, molar_composition_of_bulk), 0.0);
            new_molar_composition_of_melt = std::min(std::max(molar_composition_of_melt, molar_composition_of_bulk), 1.0);

            if (std::abs(molar_composition_of_melt - molar_composition_of_solid) > std::numeric_limits<double>::min())
              melt_molar_fraction = std::min(std::max((molar_composition_of_bulk - molar_composition_of_solid) / (molar_composition_of_melt - molar_composition_of_solid), 0.0), 1.0);
            // If solid and melt composition are the same, there is no two-phase region.
            // If we are not above the liquidus, we are below the solidus.
            else
              melt_molar_fraction = molar_volatiles_in_bulk;

            if (molar_volatiles_in_bulk > 0.0 && std::abs(molar_composition_of_bulk - molar_composition_of_solid) > std::numeric_limits<double>::min())
              molar_volatiles_in_melt = molar_volatiles_in_bulk*(molar_composition_of_melt - molar_composition_of_solid)/(molar_composition_of_bulk - molar_composition_of_solid);
          }

        return melt_molar_fraction;
      }
    }



    template <int dim>
    void
    MeltBoukare<dim>::
    melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                    std::vector<double> &melt_fractions) const
    {
      const unsigned int Fe_solid_idx = this->introspection().compositional_index_for_name("molar_Fe_in_solid");
      unsigned int Fe_melt_idx = numbers::invalid_unsigned_int;
      unsigned int porosity_idx = numbers::invalid_unsigned_int;

      if (this->include_melt_transport())
        {
          Fe_melt_idx = this->introspection().compositional_index_for_name("molar_Fe_in_melt");
          porosity_idx = this->introspection().compositional_index_for_name("porosity");
        }

      const unsigned int n_endmembers = endmember_names.size();
      EndmemberProperties endmembers(n_endmembers);

      for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
        {
          std::vector<double> endmember_mole_fractions_per_phase(n_endmembers);

          fill_endmember_properties(in, q, endmembers);

          // We need the compositions of all phases.
          double solid_composition = in.composition[q][Fe_solid_idx];
          double melt_composition = 0.0, melt_molar_fraction = 0.0;
          double bridgmanite_molar_fraction_in_solid;

          if (this->include_melt_transport())
            melt_composition = in.composition[q][Fe_melt_idx];

          convert_composition_to_fraction_of_endmembers(in.temperature[q],
                                                        solid_composition,
                                                        melt_composition,
                                                        endmembers.gibbs_energies,
                                                        endmember_mole_fractions_per_phase,
                                                        bridgmanite_molar_fraction_in_solid);


          if (this->include_melt_transport())
            melt_molar_fraction = compute_melt_molar_fraction(in.composition[q][porosity_idx],
                                                              bridgmanite_molar_fraction_in_solid,
                                                              endmembers,
                                                              endmember_mole_fractions_per_phase);

          const double solid_molar_fraction = 1.0 - melt_molar_fraction;
          const double bulk_composition = melt_composition * melt_molar_fraction + solid_composition * solid_molar_fraction;
          double molar_volatiles_in_melt = 0.0;

          const double eq_melt_molar_fraction = this->melt_fraction(in.temperature[q],
                                                                    this->get_adiabatic_conditions().pressure(in.position[q]),
                                                                    bulk_composition,
                                                                    molar_volatiles_in_melt,
                                                                    solid_composition,
                                                                    melt_composition);

          // We have to compute the endmember fractions again here because the porosity is now different.
          convert_composition_to_fraction_of_endmembers(in.temperature[q],
                                                        solid_composition,
                                                        melt_composition,
                                                        endmembers.gibbs_energies,
                                                        endmember_mole_fractions_per_phase,
                                                        bridgmanite_molar_fraction_in_solid);


          // convert from melt molar fraction to porosity
          const double solid_molar_volume = bridgmanite_molar_fraction_in_solid * endmember_mole_fractions_per_phase[febdg_idx] * endmembers.volumes[febdg_idx]
                                            + bridgmanite_molar_fraction_in_solid * endmember_mole_fractions_per_phase[mgbdg_idx] * endmembers.volumes[mgbdg_idx]
                                            + (1. - bridgmanite_molar_fraction_in_solid) * endmember_mole_fractions_per_phase[per_idx] * endmembers.volumes[per_idx]
                                            + (1. - bridgmanite_molar_fraction_in_solid) * endmember_mole_fractions_per_phase[wus_idx] * endmembers.volumes[wus_idx];
          const double melt_molar_volume = endmember_mole_fractions_per_phase[mgmelt_idx] * endmembers.volumes[mgmelt_idx]
                                           + endmember_mole_fractions_per_phase[femelt_idx] * endmembers.volumes[femelt_idx];

          melt_fractions[q] = eq_melt_molar_fraction * melt_molar_volume
                              / (eq_melt_molar_fraction * melt_molar_volume + (1.0 - eq_melt_molar_fraction) * solid_molar_volume);
        }
      return;
    }



    template <int dim>
    void
    MeltBoukare<dim>::
    evaluate(const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out) const
    {
      ReactionRateOutputs<dim> *reaction_rate_out = out.template get_additional_output<ReactionRateOutputs<dim>>();
      MeltOutputs<dim> *melt_out = out.template get_additional_output<MeltOutputs<dim>>();
      BoukareOutputs<dim> *boukare_out = out.template get_additional_output<BoukareOutputs<dim>>();
      EnthalpyOutputs<dim> *enthalpy_out = out.template get_additional_output<EnthalpyOutputs<dim>>();

      const unsigned int Fe_solid_idx = this->introspection().compositional_index_for_name("molar_Fe_in_solid");
      unsigned int Fe_melt_idx = numbers::invalid_unsigned_int;
      unsigned int porosity_idx = numbers::invalid_unsigned_int;

      if (this->include_melt_transport())
        {
          Fe_melt_idx = this->introspection().compositional_index_for_name("molar_Fe_in_melt");
          porosity_idx = this->introspection().compositional_index_for_name("porosity");
        }

      // If the temperature or pressure are zero, this model does not work.
      // This should only happen when setting the melt constraints before we have the initial temperature.
      // In this case, just fill the permeabilities and fluid viscosities and return.
      const unsigned int n_points = in.n_evaluation_points();
      for (unsigned int q=0; q<n_points; ++q)
        {
          if (in.temperature[q] == 0.0)
            {
              if (melt_out != nullptr)
                {
                  for (unsigned int q=0; q<n_points; ++q)
                    {
                      const double porosity = std::max(in.composition[q][porosity_idx],0.0);
                      melt_out->fluid_viscosities[q] = eta_f;
                      melt_out->permeabilities[q] = reference_permeability * Utilities::fixed_power<3>(porosity) * Utilities::fixed_power<2>(1.0-porosity);
                    }
                }
              return;
            }
        }

      // We need the reaction step here to conserve bulk composition.
      double reaction_time_step_size = 1.0;
      double reaction_fraction = 0.0;
      if (this->simulator_is_past_initialization())
        {
          const unsigned int number_of_reaction_steps = std::max(static_cast<unsigned int>(this->get_timestep() / this->get_parameters().reaction_time_step),
                                                                 std::max(this->get_parameters().reaction_steps_per_advection_step,1U));
          reaction_time_step_size = this->get_timestep() / static_cast<double>(number_of_reaction_steps);
          reaction_fraction = reaction_time_step_size / melting_time_scale;
        }

      const unsigned int n_endmembers = endmember_names.size();
      EndmemberProperties endmembers(n_endmembers);

      for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
        {
          std::vector<double> endmember_mole_fractions_per_phase(n_endmembers);
          std::vector<double> endmember_mole_fractions_in_composite(n_endmembers);

          fill_endmember_properties(in, q, endmembers);

          // We need the compositions of all phases.
          const double solid_composition = in.composition[q][Fe_solid_idx];
          double melt_composition = 0.0;
          double melt_molar_fraction = 0.0;
          double bridgmanite_molar_fraction_in_solid;

          if (this->include_melt_transport())
            melt_composition = in.composition[q][Fe_melt_idx];

          convert_composition_to_fraction_of_endmembers(in.temperature[q],
                                                        solid_composition,
                                                        melt_composition,
                                                        endmembers.gibbs_energies,
                                                        endmember_mole_fractions_per_phase,
                                                        bridgmanite_molar_fraction_in_solid);

          if (this->include_melt_transport())
            melt_molar_fraction = compute_melt_molar_fraction(in.composition[q][porosity_idx],
                                                              bridgmanite_molar_fraction_in_solid,
                                                              endmembers,
                                                              endmember_mole_fractions_per_phase);

          const double solid_molar_fraction = 1.0 - melt_molar_fraction;

          // Compute endmember molar fractions in the composite.
          std::vector<double> phase_mole_fractions_in_composite(n_endmembers);
          double solid_molar_mass = 0.0;
          double melt_molar_mass = 0.0;
          double melt_molar_volume = 0.0;
          double solid_molar_volume = 0.0;

          for (unsigned int i=0; i<n_endmembers; ++i)
            {
              // We assume here that the first 2 endmembers are the (solid) bridgmanite endmembers,
              // and the 3rd and 4th are the (solid) ferropericlase endmembers. The endmember_phase_fraction_in_solid
              // variable therefore indicates the fraction of bridgmanite or ferropericlase in the solid
              // (depending on which phase the endmember belongs to).
              // TODO: allow for more flexibility by letting the user set a phase vector in the input file
              // (then we can check if phase == 'bridgmanite').
              const double endmember_phase_fraction_in_solid = i<2 ? bridgmanite_molar_fraction_in_solid : (1.0 - bridgmanite_molar_fraction_in_solid);

              if (endmember_states[i] == EndmemberState::solid)
                {
                  // There are two phases in the solid (bridgmanite and ferropericlase). The mole fraction of each
                  // phase in the composite is the product of the solid fraction in the composite and the fraction
                  // of that individual phase in the solid.
                  phase_mole_fractions_in_composite[i] = solid_molar_fraction * endmember_phase_fraction_in_solid;
                  endmember_mole_fractions_in_composite[i] = phase_mole_fractions_in_composite[i] * endmember_mole_fractions_per_phase[i];
                  solid_molar_mass += endmember_mole_fractions_in_composite[i] * molar_masses[i];
                  solid_molar_volume += endmember_mole_fractions_in_composite[i] * endmembers.volumes[i];
                }
              else if (endmember_states[i] == EndmemberState::melt)
                {
                  // There is only one phase in the melt, so the mole fraction of that phase equals the melt fraction.
                  phase_mole_fractions_in_composite[i] = melt_molar_fraction;
                  endmember_mole_fractions_in_composite[i] = phase_mole_fractions_in_composite[i] * endmember_mole_fractions_per_phase[i];
                  melt_molar_mass += endmember_mole_fractions_in_composite[i] * molar_masses[i];
                  melt_molar_volume += endmember_mole_fractions_in_composite[i] * endmembers.volumes[i];
                }
              else
                AssertThrow (false, ExcNotImplemented());
            }

          const double total_molar_mass = melt_molar_mass + solid_molar_mass;
          const double total_volume = melt_molar_volume + solid_molar_volume;

          out.thermal_expansion_coefficients[q] = 0.0;
          out.specific_heat[q] = 0.0;
          out.compressibilities[q] = 0.0;

          // Average the individual endmember properties
          for (unsigned int i=0; i<n_endmembers; ++i)
            {
              out.thermal_expansion_coefficients[q] += (endmember_mole_fractions_in_composite[i] * endmembers.volumes[i] * endmembers.thermal_expansivities[i]) / total_volume;
              out.specific_heat[q] += endmember_mole_fractions_in_composite[i] * endmembers.heat_capacities[i] / total_molar_mass;


              if (endmember_states[i] == EndmemberState::solid && solid_molar_volume > 0.0)
                {
                  out.compressibilities[q] += (endmember_mole_fractions_in_composite[i] * endmembers.volumes[i])
                                              / (solid_molar_volume * endmembers.bulk_moduli[i]);
                }
            }

          if (solid_molar_volume > 0.0)
            out.densities[q] = solid_molar_mass / solid_molar_volume;
          else
            out.densities[q] = melt_molar_mass / melt_molar_volume;

          if (melt_out != nullptr)
            {
              double melt_compressiblity = 0.0;
              for (unsigned int i=0; i<n_endmembers; ++i)
                if (endmember_states[i] == EndmemberState::melt && melt_molar_volume > 0.0)
                  melt_compressiblity += (endmember_mole_fractions_in_composite[i] * endmembers.volumes[i]) / (melt_molar_volume * endmembers.bulk_moduli[i]);

              if (melt_molar_volume > 0.0)
                melt_out->fluid_densities[q] = melt_molar_mass / melt_molar_volume;
              else
                {
                  // make sure we have a useful melt density even if there is no melt to avoid density jumps
                  double mass = 0.0;
                  double volume = 0.0;
                  for (unsigned int i=0; i<n_endmembers; ++i)
                    if (endmember_states[i] == EndmemberState::melt)
                      {
                        mass += endmember_mole_fractions_per_phase[i] * molar_masses[i];
                        volume += endmember_mole_fractions_per_phase[i] * endmembers.volumes[i];
                      }
                  melt_out->fluid_densities[q] = mass/volume;
                }

              // This does not take into account the volume change due to thermal expansion of melt
              melt_out->fluid_density_gradients[q] = melt_out->fluid_densities[q] * melt_out->fluid_densities[q]
                                                     * melt_compressiblity
                                                     * this->get_gravity_model().gravity_vector(in.position[q]);
            }


          if (this->include_melt_transport())
            {
              // TODO: Alternatively, this could also be done by summing over endmember_mole_fractions_in_composite
              const double old_porosity = in.composition[q][porosity_idx];
              const double old_solid_composition = in.composition[q][Fe_solid_idx];
              const double old_melt_composition = in.composition[q][Fe_melt_idx];
              const double old_melt_molar_fraction = melt_molar_fraction;

              // in this simple model, the bulk composition is just one number, namely
              // the molar fraction of the combined iron endmembers
              const double bulk_composition = old_melt_composition * melt_molar_fraction + old_solid_composition * solid_molar_fraction;
              double molar_volatiles_in_melt = 0.0;
              double solid_composition, melt_composition;

              // calculate the melting rate as difference between the equilibrium melt fraction
              // and the solution of the previous time step, and also update melt and solid composition
              melt_molar_fraction = melt_fraction(in.temperature[q],
                                                  this->get_adiabatic_conditions().pressure(in.position[q]),
                                                  bulk_composition,
                                                  molar_volatiles_in_melt,
                                                  solid_composition,
                                                  melt_composition);

              if (boukare_out != nullptr)
                {
                  boukare_out->bulk_composition[q] = bulk_composition;
                  boukare_out->molar_volatiles_in_melt[q] = molar_volatiles_in_melt;
                }

              // We have to compute the update to the melt fraction in such a way that the bulk composition is conserved.
              const double change_of_melt_composition = reaction_fraction * assert_update_is_within_0_and_1(old_melt_composition, melt_composition - old_melt_composition);
              const double change_of_melt_fraction = reaction_fraction * assert_update_is_within_0_and_1(old_melt_molar_fraction, melt_molar_fraction - old_melt_molar_fraction);

              melt_composition = old_melt_composition + change_of_melt_composition;
              melt_molar_fraction = old_melt_molar_fraction + change_of_melt_fraction;

              if (melt_molar_fraction > 0.0 && melt_molar_fraction < 1.0)
                solid_composition = (bulk_composition - melt_molar_fraction * melt_composition)
                                    / (1.0 - melt_molar_fraction);
              const double change_of_solid_composition = solid_composition - old_solid_composition;


              // We have to compute the endmember fractions again here because the porosity is now different.
              convert_composition_to_fraction_of_endmembers(in.temperature[q],
                                                            solid_composition,
                                                            melt_composition,
                                                            endmembers.gibbs_energies,
                                                            endmember_mole_fractions_per_phase,
                                                            bridgmanite_molar_fraction_in_solid);


              // convert from melt molar fraction to porosity
              const double solid_molar_volume = bridgmanite_molar_fraction_in_solid * endmember_mole_fractions_per_phase[febdg_idx] * endmembers.volumes[febdg_idx]
                                                + bridgmanite_molar_fraction_in_solid * endmember_mole_fractions_per_phase[mgbdg_idx] * endmembers.volumes[mgbdg_idx]
                                                + (1. - bridgmanite_molar_fraction_in_solid) * endmember_mole_fractions_per_phase[per_idx] * endmembers.volumes[per_idx]
                                                + (1. - bridgmanite_molar_fraction_in_solid) * endmember_mole_fractions_per_phase[wus_idx] * endmembers.volumes[wus_idx];
              const double melt_molar_volume = endmember_mole_fractions_per_phase[mgmelt_idx] * endmembers.volumes[mgmelt_idx]
                                               + endmember_mole_fractions_per_phase[femelt_idx] * endmembers.volumes[femelt_idx];

              const double new_porosity = melt_molar_fraction * melt_molar_volume
                                          / (melt_molar_fraction * melt_molar_volume + (1.0 - melt_molar_fraction) * solid_molar_volume);
              const double porosity_change = new_porosity - old_porosity;

              // For this simple model, we only track the iron in the solid (bridgmanite) and the iron in the melt
              for (unsigned int c=0; c<in.composition[q].size(); ++c)
                {
                  // fill reaction rate outputs
                  if (reaction_rate_out != nullptr)
                    {
                      if (!include_melting_and_freezing)
                        reaction_rate_out->reaction_rates[q][c] = 0.0;
                      else if (c == Fe_solid_idx && this->get_timestep_number() > 0)
                        reaction_rate_out->reaction_rates[q][c] = change_of_solid_composition / reaction_time_step_size;
                      else if (c == Fe_melt_idx && this->get_timestep_number() > 0)
                        reaction_rate_out->reaction_rates[q][c] = change_of_melt_composition / reaction_time_step_size;
                      else if (c == porosity_idx && this->get_timestep_number() > 0)
                        reaction_rate_out->reaction_rates[q][c] = porosity_change / reaction_time_step_size;
                      else
                        reaction_rate_out->reaction_rates[q][c] = 0.0;
                    }
                  out.reaction_terms[q][c] = 0.0;
                }

              if (enthalpy_out != nullptr)
                {
                  const double melt_molar_mass = endmember_mole_fractions_per_phase[mgmelt_idx] * molar_masses[mgmelt_idx]
                                                 + endmember_mole_fractions_per_phase[femelt_idx] * molar_masses[femelt_idx];
                  const double Fe_enthalpy_of_fusion = Fe_mantle_melting_temperature * Fe_mantle_melting_entropy
                                                       + (this->get_adiabatic_conditions().pressure(in.position[q]) - melting_reference_pressure) * Fe_mantle_melting_volume;
                  const double Mg_enthalpy_of_fusion = Mg_mantle_melting_temperature * Mg_mantle_melting_entropy
                                                       + (this->get_adiabatic_conditions().pressure(in.position[q]) - melting_reference_pressure) * Mg_mantle_melting_volume;
                  double enthalpy_of_fusion = Fe_enthalpy_of_fusion * bulk_composition + Mg_enthalpy_of_fusion * (1.0-bulk_composition);
                  enthalpy_of_fusion /= melt_molar_mass;

                  enthalpy_out->enthalpies_of_fusion[q] = enthalpy_of_fusion;
                }

              // cutoff for viscosity at 30%
              const double porosity = std::min(0.3, std::max(in.composition[q][porosity_idx],0.0));
              out.viscosities[q] = (1.0 - porosity) * eta_0 * std::exp(- alpha_phi * porosity);
            }
          else
            {
              out.viscosities[q] = eta_0;

              // no melting/freezing is used in the model --> set all reactions to zero
              for (unsigned int c=0; c<in.composition[q].size(); ++c)
                {
                  if (reaction_rate_out != nullptr)
                    reaction_rate_out->reaction_rates[q][c] = 0.0;
                }
            }

          out.entropy_derivative_pressure[q]    = 0.0;
          out.entropy_derivative_temperature[q] = 0.0;
          out.thermal_conductivities[q] = thermal_conductivity;

          for (unsigned int c=0; c<in.composition[q].size(); ++c)
            out.reaction_terms[q][c] = 0.0;

          const double delta_temp = in.temperature[q]-this->get_adiabatic_conditions().temperature(in.position[q]);
          const double viscosity_bound = 1.e8;
          const double visc_temperature_dependence = std::max(std::min(std::exp(-thermal_viscosity_exponent*delta_temp/this->get_adiabatic_conditions().temperature(in.position[q])),viscosity_bound),1./viscosity_bound);
          out.viscosities[q] *= visc_temperature_dependence;
        }

      // fill melt outputs if they exist
      if (melt_out != nullptr)
        {
          const unsigned int n_points = in.n_evaluation_points();
          for (unsigned int q=0; q<n_points; ++q)
            {
              double porosity = std::max(in.composition[q][porosity_idx],0.0);

              melt_out->fluid_viscosities[q] = eta_f;
              melt_out->permeabilities[q] = reference_permeability * Utilities::fixed_power<3>(porosity) * Utilities::fixed_power<2>(1.0-porosity);

              // limit porosity to disaggregation threshold
              porosity = std::min(0.3, porosity);

              const double porosity_threshold = 0.01 * std::pow(this->get_melt_handler().melt_parameters.melt_scaling_factor_threshold, 1./3.);
              melt_out->compaction_viscosities[q] = (1.0 - porosity) * xi_0 / std::max(porosity, porosity_threshold);

              const double delta_temp = in.temperature[q]-this->get_adiabatic_conditions().temperature(in.position[q]);
              const double compaction_viscosity_bound = 1.e4;
              const double visc_temperature_dependence = std::max(std::min(std::exp(-thermal_bulk_viscosity_exponent*delta_temp/this->get_adiabatic_conditions().temperature(in.position[q])),compaction_viscosity_bound),1./compaction_viscosity_bound);

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
                             Patterns::Double (),
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
          prm.declare_entry ("Fe mantle melting temperature", "3424.5",
                             Patterns::Double(),
                             "The melting temperature of one of the components in the melting "
                             "model, the Fe mantle endmember."
                             "Units: K.");
          prm.declare_entry ("Mg mantle melting temperature", "4821.2",
                             Patterns::Double(),
                             "The melting temperature of one of the components in the melting "
                             "model, the Mg mantle endmember."
                             "Units: K.");
          prm.declare_entry ("Fe number of moles", "0.48",
                             Patterns::Double(),
                             "The number of moles of Fe atoms mixing on a pseudosite in the "
                             "mantle lattice, This is needed because we use an empirical model "
                             "fitting the full Boukare model, and can be changed to reflect "
                             "partition coefficients from other sources."
                             "Units: none.");
          prm.declare_entry ("Mg number of moles", "0.62",
                             Patterns::Double(),
                             "The number of moles of Mg atoms mixing on a pseudosite in the "
                             "mantle lattice, This is needed because we use an empirical model "
                             "fitting the full Boukare model, and can be changed to reflect "
                             "partition coefficients from other sources."
                             "Units: none.");
          prm.declare_entry ("Reference temperature", "298.15",
                             Patterns::Double(),
                             "Reference temperature used to compute the material properties"
                             "of the different endmember components."
                             "Units: K.");
          prm.declare_entry ("Reference pressure", "1e11",
                             Patterns::Double(),
                             "Reference pressure used to compute the material properties"
                             "of the different endmember components."
                             "Units: Pa.");

          prm.declare_entry ("Endmember names", "FeSiO3_bridgmanite, MgSiO3_bridgmanite, FeO_periclase, MgO_periclase, FeO_melt, MgO_melt, SiO2_melt",
                             Patterns::List(Patterns::MultipleSelection("MgSiO3_bridgmanite|FeSiO3_bridgmanite|MgO_periclase|FeO_periclase|MgO_melt|FeO_melt|SiO2_melt")),
                             "Names of the endmember components used in the equation of state and the melting model, "
                             "and whose parameters are determined by the other input parameters of this material model. "
                             "The order the parameters are given in has to be the same as the order the endmember names "
                             "are given in. "
                             "Units: none.");
          prm.declare_entry ("Endmember states", "solid, solid, solid, solid, melt, melt, melt",
                             Patterns::List(Patterns::MultipleSelection("solid|melt")),
                             "States of the endmember components used in the equation of state and the melting model. "
                             "For each endmember, this list has to define if they belong to the melt or to the solid. "
                             "The order the states are given in has to be the same as the order the 'Endmember names' "
                             "are given in. "
                             "Units: none.");

          prm.declare_entry ("Molar masses", "0.1319287, 0.1003887, 0.0718444, 0.0403044, 0.0707624708, 0.048592178, 0.048592178",
                             Patterns::List(Patterns::Double(0)),
                             "Molar masses of the different endmembers"
                             "Units: kg/mol.");
          prm.declare_entry ("Number of atoms", "5.0, 5.0, 2.0, 2.0, 2.092, 2.419, 2.419",
                             Patterns::List(Patterns::Double(0)),
                             "Number of atoms per in the formula of each endmember."
                             "Units: none.");
          prm.declare_entry ("Reference volumes", "2.534e-05, 2.445e-05, 1.206e-05, 1.125e-05, 1.2325484447664221e-05, 1.218e-05, 1.218e-05",
                             Patterns::List(Patterns::Double(0)),
                             "Reference volumes of the different endmembers."
                             "Units: $m^3$.");
          prm.declare_entry ("Reference thermal expansivities", "1.87e-05, 1.87e-05, 3.22e-05, 3.11e-05, 2.9614332469401705e-05, 2.06e-05, 2.06e-05",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal expansivities for each different endmember at the reference temperature "
                             "and reference pressure."
                             "Units: 1/K.");
          prm.declare_entry ("Reference bulk moduli", "2.81e11, 2.51e+11, 1.52e11, 1.616e11, 166652774642.11273, 2.317e11, 2.317e11",
                             Patterns::List(Patterns::Double(0)),
                             "List of bulk moduli for each different endmember at the reference temperature "
                             "and reference pressure."
                             "Units: Pa.");
          prm.declare_entry ("First derivatives of the bulk modulus", "4.14, 4.14, 4.9, 3.95, 5.0802472229003905, 4.25, 4.25",
                             Patterns::List(Patterns::Double()),
                             "The pressure derivative of the bulk modulus at the reference temperature and reference "
                             "pressure for each different endmember component."
                             "Units: none.");
          prm.declare_entry ("Second derivatives of the bulk modulus", "-1.6e-11, -1.6e-11, -3.2e-11, -2.4e-11, -3.9742163085937504e-11, -2.14e-11, -2.14e-11",
                             Patterns::List(Patterns::Double()),
                             "The second pressure derivative of the bulk modulus at the reference temperature and reference "
                             "pressure for each different endmember component."
                             "Units: 1/Pa.");
          prm.declare_entry ("Einstein temperatures", "418.1, 561.0, 297.6, 540.2, 505.75, 558.1, 558.1",
                             Patterns::List(Patterns::Double(0)),
                             "List of Einstein temperatures for each different endmember."
                             "Units: K.");
          prm.declare_entry ("Reference enthalpies", "-1082910.0, -1442310.0, -262240.0, -601570.0, -195245.49100022088, -538009.8, -538009.8",
                             Patterns::List(Patterns::Double()),
                             "List of enthalpies at the reference temperature and reference "
                             "pressure for each different endmember component."
                             "Units: J/mol.");
          prm.declare_entry ("Reference entropies", "95.0, 62.6, 58.6, 26.5, 95.0299295525918, 64.9, 64.9",
                             Patterns::List(Patterns::Double(0)),
                             "List of entropies at the reference temperature and reference "
                             "pressure for each different endmember component."
                             "Units: J/K/mol.");
          prm.declare_entry ("Reference specific heat capacities", "139.546209, 161.546581, 52.0016403, 73.1147154, 79.5326013, 79.5326013, 79.5326013",
                             Patterns::List(Patterns::Double(0)),
                             "List of specific heat capacities for each different endmember at the reference temperature "
                             "and reference pressure."
                             "Units: J/kg/K.");
          prm.declare_entry ("Linear coefficients for specific heat polynomial", "6.36191292e-03, -3.31714290e-03, 3.36163516e-03, -6.35318887e-03, -2.41909947e-03, -2.41909947e-03, -2.41909947e-03",
                             Patterns::List(Patterns::Double()),
                             "The first of three coefficients that are used to compute the specific heat capacities for each "
                             "different endmember at the reference temperature and reference pressure. "
                             "This coefficient describes the linear part of the temperature dependence. "
                             "Units: J/kg/K/K.");
          prm.declare_entry ("Second coefficients for specific heat polynomial", "-4.13886524e+06, -3.57533814e+06, -1.19540964e+06, -7.33679285e+05, -1.61692272e+06, -1.61692272e+06, -1.61692272e+06",
                             Patterns::List(Patterns::Double()),
                             "The second of three coefficients that are used to compute the specific heat capacities for each "
                             "different endmember at the reference temperature and reference pressure. This coefficient describes "
                             "the part of the temperature dependence that scales as the inverse of the square of the temperature. "
                             "Units: J K/kg.");
          prm.declare_entry ("Third coefficients for specific heat polynomial", "-464.775577, -1112.54791, 25.5067110, -592.994207, -562.222634, -562.222634, -562.222634",
                             Patterns::List(Patterns::Double()),
                             "The third of three coefficients that are used to compute the specific heat capacities for each "
                             "different endmember at the reference temperature and reference pressure. This coefficient describes "
                             "the part of the temperature dependence that scales as the inverse of the square root of the temperature"
                             "Units: J/kg/sqrt(K).");
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
          Fe_mantle_melting_temperature     = prm.get_double ("Fe mantle melting temperature");
          Mg_mantle_melting_temperature     = prm.get_double ("Mg mantle melting temperature");
          Fe_number_of_moles                = prm.get_double ("Fe number of moles");
          Mg_number_of_moles                = prm.get_double ("Mg number of moles");

          reference_temperature             = prm.get_double ("Reference temperature");
          reference_pressure                = prm.get_double ("Reference pressure");

          if (this->convert_output_to_years() == true)
            melting_time_scale *= year_in_seconds;

          AssertThrow(this->get_parameters().use_operator_splitting &&
                      this->get_parameters().reaction_solver_type == Parameters<dim>::ReactionSolverType::fixed_step,
                      ExcMessage("The melt boukare material model has to be used with operator splitting, "
                                 "and the reaction solver needs to be `fixed step'."));

          AssertThrow(melting_time_scale >= this->get_parameters().reaction_time_step,
                      ExcMessage("The reaction time step " + Utilities::to_string(this->get_parameters().reaction_time_step)
                                 + " in the operator splitting scheme is too large to compute melting rates! "
                                 "You have to choose it in such a way that it is smaller than the 'Melting time scale for "
                                 "operator splitting' chosen in the material model, which is currently "
                                 + Utilities::to_string(melting_time_scale) + "."));
          AssertThrow(melting_time_scale > 0.0,
                      ExcMessage("The Melting time scale for operator splitting must be larger than 0!"));

          // Equation of state parameters
          endmember_names = Utilities::split_string_list(prm.get("Endmember names"));
          AssertThrow(Utilities::has_unique_entries(endmember_names),
                      ExcMessage("The list of strings for the parameter "
                                 "'Material model/Melt boukare/Endmember names' contains entries more than once. "
                                 "This is not allowed. Please check your parameter file."));

          const unsigned int n_endmembers = endmember_names.size();
          AssertThrow(n_endmembers == 7 || (!this->include_melt_transport() && n_endmembers == 4),
                      ExcMessage("The list of strings for the parameter "
                                 "'Material model/Melt boukare/Endmember names' does not contain the correct "
                                 "number of entries. Please check your parameter file."));

          std::vector<std::string> endmember_states_string = Utilities::split_string_list(prm.get("Endmember states"));
          endmember_states.resize(endmember_states_string.size());

          for (unsigned int i=0; i<endmember_states_string.size(); ++i)
            {
              if (endmember_states_string[i] == "solid")
                endmember_states[i] = EndmemberState::solid;
              else if (endmember_states_string[i] == "melt")
                endmember_states[i] = EndmemberState::melt;
              else
                AssertThrow (false, ExcNotImplemented());
            }

          Utilities::MapParsing::Options options(endmember_names, "");
          options.property_name = "Molar masses";
          molar_masses = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

          options.property_name = "Number of atoms";
          number_of_atoms = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

          options.property_name = "Reference volumes";
          reference_volumes = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

          options.property_name = "Reference thermal expansivities";
          reference_thermal_expansivities = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

          options.property_name = "Reference bulk moduli";
          reference_bulk_moduli = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

          options.property_name = "First derivatives of the bulk modulus";
          bulk_modulus_pressure_derivatives = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

          options.property_name = "Second derivatives of the bulk modulus";
          bulk_modulus_second_pressure_derivatives = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

          options.property_name = "Einstein temperatures";
          Einstein_temperatures = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

          options.property_name = "Reference enthalpies";
          reference_enthalpies = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

          options.property_name = "Reference entropies";
          reference_entropies = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

          options.property_name = "Reference specific heat capacities";
          reference_specific_heats = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

          options.property_name = "Linear coefficients for specific heat polynomial";
          specific_heat_linear_coefficients = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

          options.property_name = "Second coefficients for specific heat polynomial";
          specific_heat_second_coefficients = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

          options.property_name = "Third coefficients for specific heat polynomial";
          specific_heat_third_coefficients = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

          // Check all lists have the correct length.
          AssertThrow(endmember_names.size() == endmember_states.size(),
                      ExcMessage("One of the lists that define the endmember parameters does not have the "
                                 "correct size.  Please check your parameter file."));


          // All of these are molar fractions
          // Check that all compositional fields we need exist
          AssertThrow(this->introspection().compositional_name_exists("molar_Fe_in_solid"),
                      ExcMessage("Material model melt boukare only works if there is a "
                                 "compositional field called 'molar_Fe_in_solid'."));

          if (this->include_melt_transport())
            {
              AssertThrow(this->introspection().compositional_name_exists("porosity"),
                          ExcMessage("Material model melt boukare only works if there is a "
                                     "compositional field called porosity."));
              AssertThrow(this->introspection().compositional_name_exists("molar_Fe_in_melt"),
                          ExcMessage("Material model melt boukare only works if there is a "
                                     "compositional field called 'molar_Fe_in_melt'."));
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
          && out.template get_additional_output<ReactionRateOutputs<dim>>() == nullptr)
        {
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::ReactionRateOutputs<dim>> (out.n_evaluation_points(), this->n_compositional_fields()));
        }

      if (out.template get_additional_output<BoukareOutputs<dim>>() == nullptr)
        {
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::BoukareOutputs<dim>> (out.n_evaluation_points()));
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
                                   "A material model that implements a simplified version of the melting "
                                   "model of Boukare et al. (https://doi.org/10.1002/2015JB011929) for the "
                                   "lowermost mantle and uses it to compute the material parameters "
                                   "required for the modeling of melt transport, including melting and "
                                   "solidification and the corresponding changes in composition."
                                   "The model parameterizes the composition (which includes the components "
                                   "MgO, FeO and SiO2) as a mixture between two endmembers (one iron-bearing "
                                   "and one magnesium-bearing). The equation of state considers three phases: "
                                   "bridgmanite, ferropericlase, and melt (each with their individual "
                                   "compositions). "
                                   "More details can be found in Dannberg, J., Myhill, R., Gassmöller, R., "
                                   "and Cottaar, S. (2021). The morphology, evolution and seismic visibility "
                                   "of partial melt at the core–mantle boundary: implications for ULVZs. "
                                   "Geophysical Journal International, 227(2), 1028-1059.")
  }
}
