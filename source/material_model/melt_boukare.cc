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
#include <aspect/gravity_model/interface.h>

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
    limit_update_to_0_and_1 (const double old_value,
                             const double change_of_value) const
    {
      if (old_value + change_of_value < 0)
        return -old_value;
      else if (old_value + change_of_value > 1)
        return 1.0 - old_value;
      else
        return change_of_value;
    }


    template <int dim>
    double
    MeltBoukare<dim>::
    melt_fraction (const double temperature,
                   const double pressure,
                   const double bulk_composition,
                   double &new_solid_composition,
                   double &new_melt_composition) const
    {
      if (temperature == 0)
        return 0;

      // TODO: fix
      return 0;

      // after Phipps Morgan, Jason. "Thermodynamics of pressure release melting of a veined plum pudding mantle."
      // Geochemistry, Geophysics, Geosystems 2.4 (2001). Values below taken from table A1.

      const double P = pressure;                   // pressure in Pa
      const double T = temperature;                // temperature in K
      const double R = constants::gas_constant;    // Ideal Gas Constant

      const double T_Mg_pyrolite_melting_surface = 2163.0;  // Kelvin at 1 atmosphere - reference melting temperature for Mg pyrolite endmember
      const double T_Fe_pyrolite_melting_surface = 1478.0;  // Kelvin at 1 atmosphere - reference melting temperature for Fe pyrolite endmember

      const double dS_Mg_pyrolite = 60.0;                   // entropy change of melting in J/mol K = 901.961 J/(kg K)
      const double dS_Fe_pyrolite = 60.0;                   // entropy change of melting in J/mol K = 623.025 J/(kg K)

      const double vliq_ref = 49.0/1.e6;                    // reference molar volume of the melt in m3/mol
      const double vsol_ref_Mg_pyrolite = 46.0/1.e6;        // reference molar volume of solid Mg pyrolite endmember in m3/mol
      const double vsol_ref_Fe_pyrolite = 46.0/1.e6;        // reference molar volume of solid Fe pyrolite endmember in m3/mol

      const double surface_pressure = 1.e5;                 // 1 atmosphere = 1e5 Pa

      // Free Energy Change Delta_G due to Melting as a function of temperature and pressure, for Mg pyrolite endmember and Fe pyrolite endmember.
      // Equation (A9) in Phipps Morgan (2001).
      const double dG_Mg_pyrolite = 0.0;
      const double dG_Fe_pyrolite = 0.0; // liquid - solid

      // Mole Fraction of Each Component in Coexisting Liquids and Solids, equations (A10 - A12) in Phipps Morgan (2001)
      double Xl_Mg_pyrolite = 1.0 - (1.0 - exp(dG_Mg_pyrolite/(2.0*R*T)))/(exp(dG_Fe_pyrolite/(2.0*R*T))-exp(dG_Mg_pyrolite/(2.0*R*T)));
      double Xs_Mg_pyrolite = Xl_Mg_pyrolite * exp(dG_Mg_pyrolite/(2.0*R*T));

      double melt_fraction;
      if (Xs_Mg_pyrolite <= bulk_composition)      // below the solidus
        {
          melt_fraction = 0;
          new_solid_composition = bulk_composition;
        }
      else if (Xl_Mg_pyrolite >= bulk_composition) // above the liquidus
        {
          melt_fraction = 1;
          new_melt_composition = bulk_composition;
        }
      else                                         // between solidus and liquidis
        {
          new_solid_composition = Xs_Mg_pyrolite;
          new_melt_composition = Xl_Mg_pyrolite;
          melt_fraction = (bulk_composition - Xs_Mg_pyrolite)/(Xl_Mg_pyrolite - Xs_Mg_pyrolite);
        }

      return melt_fraction;
    }


    template <int dim>
    void
    MeltBoukare<dim>::
    melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                    std::vector<double> &melt_fractions) const
    {
      double bulk_composition = 0.89;
      double solid_composition = 0.89;
      double melt_composition = 0.89;

      for (unsigned int q=0; q<in.temperature.size(); ++q)
        {
          // TODO
          melt_fractions[q] = this->melt_fraction(in.temperature[q],
                                                  std::max(0.0, in.pressure[q]),
                                                  bulk_composition,
                                                  solid_composition,
                                                  melt_composition);
        }
      return;
    }


    template <int dim>
    void
    MeltBoukare<dim>::
    evaluate(const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out) const
    {
      ReactionRateOutputs<dim> *reaction_rate_out = out.template get_additional_output<ReactionRateOutputs<dim> >();
      MeltOutputs<dim> *melt_out = out.template get_additional_output<MeltOutputs<dim> >();

      // if the temperature or pressure are zero, this model does not work
      // this should only happen when wetting the melt constraints before we have the initial temperature
      // in this case, just fill the permeabilities and fluid viscosities and return
      for (unsigned int q=0; q<in.position.size(); ++q)
      {
        if(in.temperature[q] == 0.0)
        {
		  if (melt_out != nullptr)
			{
			  const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
			  for (unsigned int q=0; q<in.position.size(); ++q)
				{
				  double porosity = std::max(in.composition[q][porosity_idx],0.0);
				  melt_out->fluid_viscosities[q] = eta_f;
				  melt_out->permeabilities[q] = reference_permeability * std::pow(porosity,3) * std::pow(1.0-porosity,2);
				}
			}
		  return;
        }
      }



      // get the indices for the different components and phases we track
      const unsigned int bridgmanite_idx = this->introspection().compositional_index_for_name("bridgmanite");
      const unsigned int Fe_in_bridgmanite_idx = this->introspection().compositional_index_for_name("iron_fraction_in_bridgmanite");
      const unsigned int Fe_in_periclase_idx = this->introspection().compositional_index_for_name("iron_fraction_in_periclase");

      // get indices of the different endmembers
      // TODO: get from the list of endmembers, for now: MgSiO3_bridgmanite, FeSiO3_bridgmanite, MgO_periclase, FeO_periclase
      // TODO: also get the opposite (name of index)
      const unsigned int febdg_idx = 0;
      const unsigned int mgbdg_idx = 1;
      const unsigned int wus_idx = 2;
      const unsigned int per_idx = 3;
      const unsigned int mgmelt_idx = 4;
      const unsigned int femelt_idx = 5;
      const unsigned int simelt_idx = 6;


      for (unsigned int q=0; q<in.position.size(); ++q)
        {
          // TODO: get a vector of endmember phases?
          const unsigned int n_endmembers = endmember_names.size();

          std::vector<double> endmember_mole_fractions_per_phase(n_endmembers);
          std::vector<double> endmember_mole_fractions_in_composite(n_endmembers);
          std::vector<double> phase_mass_fractions(n_endmembers);

          double porosity = 0.0;

          // the two solid phases are bridgmanite and periclase
          // all mass fractions represent the fraction of mass in the solid phase
          const double bridgmanite_mass_fraction = in.composition[q][bridgmanite_idx];
          phase_mass_fractions[mgbdg_idx] = bridgmanite_mass_fraction;
          phase_mass_fractions[febdg_idx] = bridgmanite_mass_fraction;

          const double periclase_mass_fraction = 1.0 - bridgmanite_mass_fraction;
          phase_mass_fractions[per_idx] = periclase_mass_fraction;
          phase_mass_fractions[wus_idx] = periclase_mass_fraction;

          // we have to track the bulk composition of both solid phases

          // TODO: remove this part
          const double molar_FeSiO3_in_bridgmanite = in.composition[q][Fe_in_bridgmanite_idx];
          const double molar_MgSiO3_in_bridgmanite = 1.0 - molar_FeSiO3_in_bridgmanite;
          const double molar_FeO_in_periclase = in.composition[q][Fe_in_periclase_idx];
          const double molar_MgO_in_periclase = 1.0 - molar_FeO_in_periclase;

          endmember_mole_fractions_per_phase[febdg_idx] = in.composition[q][Fe_in_bridgmanite_idx];
          endmember_mole_fractions_per_phase[mgbdg_idx] = 1.0 - endmember_mole_fractions_per_phase[febdg_idx];
          endmember_mole_fractions_per_phase[wus_idx] = in.composition[q][Fe_in_periclase_idx];
          endmember_mole_fractions_per_phase[per_idx] = 1.0 - endmember_mole_fractions_per_phase[wus_idx];

          if (this->include_melt_transport())
            {
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
              const unsigned int Mg_in_melt_idx = this->introspection().compositional_index_for_name("magnesium_fraction_in_the_melt");
              const unsigned int Fe_in_melt_idx = this->introspection().compositional_index_for_name("iron_fraction_in_the_melt");

              porosity = in.composition[q][porosity_idx];
              endmember_mole_fractions_per_phase[mgmelt_idx] = in.composition[q][Mg_in_melt_idx];
              endmember_mole_fractions_per_phase[femelt_idx] = in.composition[q][Fe_in_melt_idx];
              endmember_mole_fractions_per_phase[simelt_idx] = 1.0 - endmember_mole_fractions_per_phase[mgmelt_idx] - endmember_mole_fractions_per_phase[femelt_idx];
            }

          out.thermal_expansion_coefficients[q] = 0.0;
          out.specific_heat[q] = 0.0;
          out.compressibilities[q] = 0.0;


          // TODO: more descriptive variable names and remove excess brackets once we have a unit test

          std::vector<double> endmember_volumes(n_endmembers);
          std::vector<double> endmember_densities(n_endmembers);
          std::vector<double> endmember_gibbs_energies(n_endmembers);
          std::vector<double> endmember_entropies(n_endmembers);
          std::vector<double> endmember_thermal_expansivities(n_endmembers);
          std::vector<double> endmember_bulk_moduli(n_endmembers);
          std::vector<double> endmember_heat_capacities(n_endmembers);

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
              double intVdP = 0.;
              double dintVdpdT = 0.;

              if (in.pressure[q] != reference_pressure && in.pressure[q] > 0.0)
                {
                  intVdP = in.pressure[q] * reference_volumes[i]
                           * (1. - a + (a * (std::pow(1. - b * Pth, 1. - c) - std::pow(1. + b * (in.pressure[q] - Pth), 1. - c))
                                        / (b * (c - 1.) * in.pressure[q])));
                  intVdP -= reference_pressure * reference_volumes[i]
                            * (1. - a + (a * (std::pow(1. - b * Pth, 1. - c) - std::pow(1. + b * (reference_pressure - Pth), 1. - c))
                                         / (b * (c - 1.) * reference_pressure)));

                  const double prefactor = reference_volumes[i] * reference_thermal_expansivities[i] * reference_bulk_moduli[i] * a * heat_capacity_ratio;
                  dintVdpdT = prefactor * (std::pow(1. + b * (in.pressure[q] - Pth), -c) - std::pow(1. + b * (reference_pressure - Pth), -c));
                }


              endmember_gibbs_energies[i] = G_Pref_Tf + intVdP;
              endmember_entropies[i] =  reference_entropies[i] + endmember_entropy_thermal_addition(in.temperature[q], i) + dintVdpdT;
              endmember_volumes[i] = reference_volumes[i]*(1 - a * (1. - std::pow(1. + b * (in.pressure[q] - Pth), -c)));
              endmember_densities[i] = molar_masses[i]/endmember_volumes[i];

              endmember_bulk_moduli[i] = reference_bulk_moduli[i] * (1. + b * (in.pressure[q] - Pth))
                                         * (a + (1. - a) * std::pow(1. + b * (in.pressure[q] - Pth), c));

              const double C_V0 = endmember_molar_heat_capacity(reference_temperature, i);
              const double C_V = endmember_molar_heat_capacity(in.temperature[q], i);
              endmember_thermal_expansivities[i] = reference_thermal_expansivities[i] * (C_V / C_V0) *
                                                   1. / ((1. + b * (in.pressure[q] - Pth)) *
                                                         (a + (1. - a) * std::pow(1 + b * (in.pressure[q] - Pth), c)));

              const double Cp_ref = reference_specific_heats[i] + specific_heat_linear_coefficients[i] * in.temperature[q]
                                    + specific_heat_second_coefficients[i] * std::pow(in.temperature[q], -2.)
                                    + specific_heat_third_coefficients[i] * std::pow(in.temperature[q], -0.5);

              const double pressure = in.pressure[q] > 0.0 ? reference_pressure : 0.0;
              const double dSdT0 = reference_volumes[i] * reference_bulk_moduli[i] * std::pow(heat_capacity_ratio * reference_thermal_expansivities[i], 2.0)
                                   * (std::pow(1. + b * (in.pressure[q] - Pth), -1.-c) - std::pow(1. + b * (pressure - Pth), -1.- c));

              const double relative_T = Einstein_temperatures[i] / in.temperature[q];
              const double dSdT = dSdT0 + dintVdpdT * (1 - 2./relative_T + 2./(std::exp(relative_T) - 1.)) * relative_T/in.temperature[q];

              endmember_heat_capacities[i] = Cp_ref + in.temperature[q] * dSdT;
            }

          // Average the individual endmember properties
          double melt_molar_volume = 0.0;
          for (unsigned int i=0; i<n_endmembers; ++i)
              if (endmember_states[i] == EndmemberState::melt)
                melt_molar_volume += endmember_mole_fractions_per_phase[i] * endmember_volumes[i];

          // convert to mole fractions
          const double n_moles_in_bridgmanite = bridgmanite_mass_fraction / (molar_MgSiO3_in_bridgmanite * molar_masses[mgbdg_idx]
                                                                             + molar_FeSiO3_in_bridgmanite * molar_masses[febdg_idx]);
          const double n_moles_in_ferropericlase = periclase_mass_fraction / (molar_MgO_in_periclase * molar_masses[per_idx]
                                                                              + molar_FeO_in_periclase * molar_masses[wus_idx]);
          const double n_moles_in_the_solid = n_moles_in_bridgmanite + n_moles_in_ferropericlase;

          const double bridgmanite_molar_fraction_in_solid = n_moles_in_bridgmanite / n_moles_in_the_solid;
          const double ferropericlase_molar_fraction_in_solid = 1.0 - bridgmanite_molar_fraction_in_solid;

          double solid_molar_volume = bridgmanite_molar_fraction_in_solid * molar_MgSiO3_in_bridgmanite * endmember_volumes[mgbdg_idx]
                                      + bridgmanite_molar_fraction_in_solid * molar_FeSiO3_in_bridgmanite * endmember_volumes[febdg_idx]
                                      + ferropericlase_molar_fraction_in_solid * molar_MgO_in_periclase * endmember_volumes[per_idx]
                                      + ferropericlase_molar_fraction_in_solid * molar_FeO_in_periclase * endmember_volumes[wus_idx];

          const double n_moles_total = melt_molar_volume > 0
                                       ?
                                       porosity / melt_molar_volume + (1.0 - porosity) / solid_molar_volume
                                       :
                                       0.0;
          double melt_molar_fraction = melt_molar_volume > 0
                                       ?
                                       porosity / (melt_molar_volume * n_moles_total)
                                       :
                                       0.0;
          const double solid_molar_fraction = 1.0 - melt_molar_fraction;

          std::vector<double> phase_mole_fractions_in_composite(n_endmembers);
          double solid_molar_mass = 0.0;
          double melt_molar_mass = 0.0;
          double total_volume = 0.0;
          // TODO: fix this using a phase vector (phase == 'bridgmanite')
          for (unsigned int i=0; i<n_endmembers; ++i)
            {
              if (i<2)
                {
                  phase_mole_fractions_in_composite[i] = solid_molar_fraction * bridgmanite_molar_fraction_in_solid;
                  solid_molar_mass += bridgmanite_molar_fraction_in_solid * endmember_mole_fractions_per_phase[i] * molar_masses[i];
                }
              else if (i>=2 && i<4)
                {
                  phase_mole_fractions_in_composite[i] = solid_molar_fraction * ferropericlase_molar_fraction_in_solid;
                  solid_molar_mass += ferropericlase_molar_fraction_in_solid * endmember_mole_fractions_per_phase[i] * molar_masses[i];
                }
              else
                {
                  phase_mole_fractions_in_composite[i] = melt_molar_fraction;
                  melt_molar_mass += endmember_mole_fractions_per_phase[i] * molar_masses[i];
                }

              endmember_mole_fractions_in_composite[i] = phase_mole_fractions_in_composite[i] * endmember_mole_fractions_per_phase[i];
              total_volume += endmember_mole_fractions_in_composite[i] * endmember_volumes[i];
            }

          for (unsigned int i=0; i<n_endmembers; ++i)
            {
              out.thermal_expansion_coefficients[q] += (endmember_mole_fractions_in_composite[i] * endmember_volumes[i] * endmember_thermal_expansivities[i]) / total_volume;
              out.specific_heat[q] += endmember_mole_fractions_in_composite[i] * endmember_heat_capacities[i];

              if (endmember_states[i] == EndmemberState::solid)
              out.compressibilities[q] += (endmember_mole_fractions_in_composite[i] * endmember_volumes[i]) / (solid_molar_volume * endmember_bulk_moduli[i]);
            }

          out.densities[q] = solid_molar_mass / solid_molar_volume;

          if (melt_out != nullptr)
            {
        	  double melt_compressiblity = 0;
              for (unsigned int i=0; i<n_endmembers; ++i)
                if (endmember_states[i] == EndmemberState::melt)
            	  melt_compressiblity += (endmember_mole_fractions_in_composite[i] * endmember_volumes[i]) / (melt_molar_volume * endmember_bulk_moduli[i]);

              melt_out->fluid_densities[q] = melt_molar_mass / melt_molar_volume;
              // TODO: this does not take into account the volume change due to thermal expansion of melt
              melt_out->fluid_density_gradients[q] = melt_out->fluid_densities[q] * melt_out->fluid_densities[q]
													 * melt_compressiblity
													 * this->get_gravity_model().gravity_vector(in.position[q]);
            }


          if (this->include_melt_transport())
            {
              // Calculate the melting rate
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

              const double old_porosity = in.composition[q][porosity_idx];
              const double old_solid_composition = endmember_mole_fractions_in_composite[febdg_idx] + endmember_mole_fractions_in_composite[wus_idx];
              const double old_melt_composition = endmember_mole_fractions_in_composite[femelt_idx];

              double solid_composition = std::min(1.0, std::max(old_solid_composition,0.0));
              double melt_composition = std::min(1.0, std::max(old_melt_composition,0.0));

              // in this simple model, the bulk composition is just one number, namely
              // the molar fraction of the combined iron endmembers
              const double bulk_composition = endmember_mole_fractions_in_composite[febdg_idx]
                                              + endmember_mole_fractions_in_composite[wus_idx]
                                              + endmember_mole_fractions_in_composite[femelt_idx];

              // calculate the melting rate as difference between the equilibrium melt fraction
              // and the solution of the previous time step, and also update melt and solid composition
              melt_molar_fraction = melt_fraction(in.temperature[q],
                                                  this->get_adiabatic_conditions().pressure(in.position[q]),
                                                  bulk_composition,
                                                  solid_composition,
                                                  melt_composition);

              // These outputs are the molar compositions
              // TODO: This is only valid if there is only one solid and one liquid phase, each with two endmembers
              solid_molar_volume = (1.0 - solid_composition) * endmember_volumes[mgbdg_idx]
                                   + solid_composition * endmember_volumes[febdg_idx];
              melt_molar_volume = (1.0 - melt_composition) * endmember_volumes[mgmelt_idx]
                                  + melt_composition * endmember_volumes[femelt_idx];

              const double new_porosity = melt_molar_fraction * melt_molar_volume
                                          / (melt_molar_fraction * melt_molar_volume + (1.0 - melt_molar_fraction) * solid_molar_volume);


              // do not allow negative porosity or porosity > 1
              const double porosity_change = limit_update_to_0_and_1(old_porosity, new_porosity - old_porosity);
              const double change_of_solid_composition = limit_update_to_0_and_1(old_solid_composition, solid_composition - old_solid_composition);
              const double change_of_melt_composition = limit_update_to_0_and_1(old_melt_composition, melt_composition - old_melt_composition);


              // For this simple model, we only track the iron in the solid (bridgmanite) and the iron in the melt
              // TODO: make this more general to work with the full equation of state
              const unsigned int solid_composition_idx = this->introspection().compositional_index_for_name("iron_fraction_in_bridgmanite");
              const unsigned int melt_composition_idx = this->introspection().compositional_index_for_name("iron_fraction_in_the_melt");

              for (unsigned int c=0; c<in.composition[q].size(); ++c)
                {
                  // fill reaction rate outputs
                  if (reaction_rate_out != NULL)
                    {
                      if (!include_melting_and_freezing)
                        reaction_rate_out->reaction_rates[q][c] = 0.0;
                      if (c == solid_composition_idx && this->get_timestep_number() > 0)
                        reaction_rate_out->reaction_rates[q][c] = change_of_solid_composition / melting_time_scale;
                      else if (c == melt_composition_idx && this->get_timestep_number() > 0)
                        reaction_rate_out->reaction_rates[q][c] = change_of_melt_composition / melting_time_scale;
                      else if (c == porosity_idx && this->get_timestep_number() > 0)
                        reaction_rate_out->reaction_rates[q][c] = porosity_change / melting_time_scale;
                      else
                        reaction_rate_out->reaction_rates[q][c] = 0.0;
                    }
                  out.reaction_terms[q][c] = 0.0;
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
                  if (reaction_rate_out != nullptr)
                    reaction_rate_out->reaction_rates[q][c] = 0.0;
                }
            }

          out.entropy_derivative_pressure[q]    = 0.0;
          out.entropy_derivative_temperature[q] = 0.0;
          out.thermal_conductivities[q] = thermal_conductivity;

          for (unsigned int c=0; c<in.composition[q].size(); ++c)
            out.reaction_terms[q][c] = 0.0;

          double visc_temperature_dependence = 1.0;
          const double delta_temp = in.temperature[q]-this->get_adiabatic_conditions().temperature(in.position[q]);
          visc_temperature_dependence = std::max(std::min(std::exp(-thermal_viscosity_exponent*delta_temp/this->get_adiabatic_conditions().temperature(in.position[q])),1e4),1e-4);
          out.viscosities[q] *= visc_temperature_dependence;
        }

      // fill melt outputs if they exist

      if (melt_out != nullptr)
        {
          const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

          for (unsigned int q=0; q<in.position.size(); ++q)
            {
              double porosity = std::max(in.composition[q][porosity_idx],0.0);

              melt_out->fluid_viscosities[q] = eta_f;
              melt_out->permeabilities[q] = reference_permeability * std::pow(porosity,3) * std::pow(1.0-porosity,2);
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

          prm.declare_entry ("Endmember names", "FeSiO3_bridgmanite, MgSiO3_bridgmanite, FeO_periclase, MgO_periclase, MgO_melt, FeO_melt, SiO2_melt",
                             Patterns::List(Patterns::MultipleSelection("MgSiO3_bridgmanite|FeSiO3_bridgmanite|MgO_periclase|FeO_periclase|MgO_melt|FeO_melt|SiO2_melt")),
                             "Names of the endmember components used in the equation of state and the melting model, "
                             "and whose parameters are determined by the other input parameters of this material model. "
                             "The order the parameters are given in has to be the same as the order the endmember names "
                             "are given in. "
                             "Units: none.");
          prm.declare_entry ("Endmember states", "solid, solid, solid, solid, melt, melt, melt",
                             Patterns::List(Patterns::MultipleSelection("solid|melt")),
                             "States of the endmember components used in the equation of state and the melting model. "
                             "For each endmember, this leist has to define if they belong to the melt or to the solid. "
                             "The order the states are given in has to be the same as the order the 'Endmember names' "
                             "are given in. "
                             "Units: none.");

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
          const unsigned int n_endmembers = this->include_melt_transport() ? 7 : 4;
          endmember_names = Utilities::split_string_list(prm.get("Endmember names"));
          AssertThrow(Utilities::has_unique_entries(endmember_names),
                      ExcMessage("The list of strings for the parameter "
                                 "'Material model/Melt boukare/Endmember names' contains entries more than once. "
                                 "This is not allowed. Please check your parameter file."));

          AssertThrow(endmember_names.size() == n_endmembers,
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
          AssertThrow(endmember_names.size() == endmember_states.size(),
                      ExcMessage("One of the lists that define the endmember parameters does not have the "
                    		  "correct size.  Please check your parameter file."));


          // All of these are molar fractions
          // Check that all compositional fields we need exist
          AssertThrow(this->introspection().compositional_name_exists("bridgmanite"),
                      ExcMessage("Material model melt boukare only works if there is a "
                                 "compositional field called bridgmanite."));
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
