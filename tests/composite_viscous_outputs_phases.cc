/*
  Copyright (C) 2024 by the authors of the ASPECT code.

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

#include <aspect/simulator.h>
#include <aspect/material_model/utilities.h>
#include <aspect/material_model/rheology/composite_visco_plastic.h>
#include <aspect/simulator_signals.h>

namespace aspect
{
  template <int dim>
  void f(const aspect::SimulatorAccess<dim> &simulator_access,
         aspect::Assemblers::Manager<dim> &)
  {
    // This function tests whether the composite creep rheology is producing
    // the correct composite viscosity and partial strain rates corresponding to
    // the different creep mechanisms incorporated into the rheology.
    // It is assumed that each individual creep mechanism has already been tested.

    using namespace aspect::MaterialModel;

    // First, we set up a few objects which are used by the rheology model.
    aspect::ParameterHandler prm;
    const std::vector<std::string> list_of_composition_names = simulator_access.introspection().get_composition_names();
    MaterialUtilities::PhaseFunction<dim> phase_function;
    phase_function.initialize_simulator (simulator_access.get_simulator());
    phase_function.declare_parameters (prm);
    prm.set("Define transition by depth instead of pressure", "false");
    prm.set("Phase transition pressures", "3e9");
    prm.set("Phase transition pressure widths", "1e9");
    prm.set("Phase transition temperatures", "273");
    prm.set("Phase transition Clapeyron slopes", "0");
    phase_function.parse_parameters (prm);

    std::vector<unsigned int> n_phases_for_each_composition = phase_function.n_phases_for_each_composition();
    // Currently, phase_function.n_phases_for_each_composition() returns a list of length
    // equal to the total number of compositions, whether or not they are chemical compositions.
    // The equation_of_state (multicomponent incompressible) requires a list only for
    // chemical compositions.
    std::vector<unsigned int> n_phases_for_each_chemical_composition = {n_phases_for_each_composition[0]};
    std::vector<unsigned int> n_phase_transitions_for_each_chemical_composition = {n_phases_for_each_composition[0] - 1};
    unsigned int n_phases = n_phases_for_each_composition[0];
    for (auto i : simulator_access.introspection().chemical_composition_field_indices())
      {
        n_phases_for_each_chemical_composition.push_back(n_phases_for_each_composition[i+1]);
        n_phase_transitions_for_each_chemical_composition.push_back(n_phases_for_each_composition[i+1] - 1);
        n_phases += n_phases_for_each_composition[i+1];
      }

    const unsigned int composition = 0;
    const std::vector<double> volume_fractions = {1.};
    std::vector<double> phase_function_values = {0.};
    const std::vector<unsigned int> n_phase_transitions_per_composition = n_phase_transitions_for_each_chemical_composition;

    // Next, we initialise instances of the composite rheology and
    // individual creep mechanisms.
    std::unique_ptr<Rheology::CompositeViscoPlastic<dim>> composite_creep;
    composite_creep = std::make_unique<Rheology::CompositeViscoPlastic<dim>>();
    composite_creep->initialize_simulator (simulator_access.get_simulator());
    composite_creep->declare_parameters(prm);
    MaterialUtilities::PhaseFunction<dim>::declare_parameters(prm);
    prm.set("Viscosity averaging scheme", "isostrain");
    prm.set("Include diffusion creep in composite rheology", "true");
    prm.set("Include dislocation creep in composite rheology", "true");
    prm.set("Include Peierls creep in composite rheology", "true");
    prm.set("Include Drucker Prager plasticity in composite rheology", "true");
    prm.set("Peierls creep flow law", "viscosity approximation");
    prm.set("Cohesions", "background:1e9|5e8");
    prm.set("Maximum yield stress", "5e10");
    composite_creep->parse_parameters(prm, std::make_unique<std::vector<unsigned int>>(n_phases_for_each_chemical_composition));

    std::unique_ptr<Rheology::DiffusionCreep<dim>> diffusion_creep;
    diffusion_creep = std::make_unique<Rheology::DiffusionCreep<dim>>();
    diffusion_creep->initialize_simulator (simulator_access.get_simulator());
    diffusion_creep->declare_parameters(prm);
    diffusion_creep->parse_parameters(prm, std::make_unique<std::vector<unsigned int>>(n_phases_for_each_chemical_composition));

    std::unique_ptr<Rheology::DislocationCreep<dim>> dislocation_creep;
    dislocation_creep = std::make_unique<Rheology::DislocationCreep<dim>>();
    dislocation_creep->initialize_simulator (simulator_access.get_simulator());
    dislocation_creep->declare_parameters(prm);
    dislocation_creep->parse_parameters(prm, std::make_unique<std::vector<unsigned int>>(n_phases_for_each_chemical_composition));

    std::unique_ptr<Rheology::PeierlsCreep<dim>> peierls_creep;
    peierls_creep = std::make_unique<Rheology::PeierlsCreep<dim>>();
    peierls_creep->initialize_simulator (simulator_access.get_simulator());
    peierls_creep->declare_parameters(prm);
    peierls_creep->parse_parameters(prm, std::make_unique<std::vector<unsigned int>>(n_phases_for_each_chemical_composition));

    std::unique_ptr<Rheology::DruckerPragerPower<dim>> drucker_prager_power;
    drucker_prager_power = std::make_unique<Rheology::DruckerPragerPower<dim>>();
    drucker_prager_power->initialize_simulator (simulator_access.get_simulator());
    drucker_prager_power->declare_parameters(prm);
    prm.set("Cohesions", "background:1e9|5e8");
    prm.set("Maximum yield stress", "5e10");
    drucker_prager_power->parse_parameters(prm, std::make_unique<std::vector<unsigned int>>(n_phases_for_each_chemical_composition));

    // The creep components are arranged in series with each other.
    // This package of components is then arranged in parallel with
    // a strain rate limiter with a constant viscosity lim_visc.
    // The whole system is then arranged in series with a viscosity limiter with
    // viscosity max_visc.
    // lim_visc is equal to (min_visc*max_visc)/(max_visc - min_visc)
    double min_visc = prm.get_double("Minimum viscosity");
    double max_visc = prm.get_double("Maximum viscosity");
    double lim_visc = (min_visc*max_visc)/(max_visc - min_visc);

    // Assign values to the variables which will be passed to compute_viscosity
    // The test involves pure shear calculations at variable pressure and temperature
    double temperature;
    double pressure;
    const double grain_size = 1.e-3;
    SymmetricTensor<2,dim> strain_rate;
    strain_rate[0][0] = -1e-11;
    strain_rate[0][1] = 0.;
    strain_rate[1][1] = 1e-11;
    strain_rate[2][0] = 0.;
    strain_rate[2][1] = 0.;
    strain_rate[2][2] = 0.;

    std::cout << "temperature (K)   phase transition progress   eta (Pas)   creep stress (Pa)   edot_ii (/s)   edot_ii fractions (diff, disl, prls, drpr, max)" << std::endl;

    // Loop through strain rates, tracking whether there is a discrepancy in
    // the decomposed strain rates.
    bool error = false;
    double viscosity;
    double total_strain_rate;
    double creep_strain_rate;
    double creep_stress;
    double diff_stress;
    double disl_stress;
    double prls_stress;
    double drpr_stress;
    std::vector<double> partial_strain_rates(5, 0.);

    for (unsigned int i=0; i <= 2; i++)
      {
        pressure = 1.e9 + i*2.e9;
        std::cout << "pressure: " << pressure / 1.e9 << " GPa" << std::endl;
        for (unsigned int j = 0; j <= 10; j++)
          {
            temperature = 1000. + j*100.;

            // Compute the phase function values
            // The depth and gravity are set to zero because they are unused
            // when phase functions are calculated by pressure.
            // The phase index is set to invalid_unsigned_int, because it is only used internally
            // in phase_average_equation_of_state_outputs to loop over all existing phases
            MaterialUtilities::PhaseFunctionInputs<dim> phase_inputs(temperature,
                                                                     pressure,
                                                                     0., 0.,
                                                                     numbers::invalid_unsigned_int);

            // Compute value of phase functions
            for (unsigned int j=0; j < phase_function.n_phase_transitions(); ++j)
              {
                phase_inputs.phase_transition_index = j;
                phase_function_values[j] = phase_function.compute_value(phase_inputs);
              }

            Rheology::DruckerPragerParameters p = drucker_prager_power->compute_drucker_prager_parameters(composition, phase_function_values, n_phase_transitions_per_composition);

            // Compute the viscosity
            viscosity = composite_creep->compute_viscosity(pressure, temperature, grain_size, volume_fractions, strain_rate, partial_strain_rates, phase_function_values, n_phase_transitions_per_composition);
            total_strain_rate = std::accumulate(partial_strain_rates.begin(), partial_strain_rates.end(), 0.);

            // The creep strain rate is calculated by subtracting the strain rate
            // of the max viscosity dashpot from the total strain rate
            // The creep stress is then calculated by subtracting the stress running
            // through the strain rate limiter from the total stress
            creep_strain_rate = total_strain_rate - partial_strain_rates[4];
            creep_stress = 2.*(viscosity*total_strain_rate - lim_visc*creep_strain_rate);

            // Print the output
            std::cout << temperature << ' ' << phase_function_values[0] << ' ' << viscosity << ' ' << creep_stress << ' ' << total_strain_rate;
            for (unsigned int i=0; i < partial_strain_rates.size(); ++i)
              {
                std::cout << ' ' << partial_strain_rates[i]/total_strain_rate;
              }
            std::cout << std::endl;

            // The following lines test that each individual creep mechanism
            // experiences the same creep stress

            // Each creep mechanism should experience the same stress
            diff_stress = 2.*partial_strain_rates[0]*diffusion_creep->compute_viscosity(pressure, temperature, grain_size, composition);
            disl_stress = 2.*partial_strain_rates[1]*dislocation_creep->compute_viscosity(partial_strain_rates[1], pressure, temperature, composition);
            prls_stress = 2.*partial_strain_rates[2]*peierls_creep->compute_viscosity(partial_strain_rates[2], pressure, temperature, composition);
            if (partial_strain_rates[3] > 0.)
              {
                drpr_stress = 2.*partial_strain_rates[3]*drucker_prager_power->compute_viscosity(p.cohesion,
                              p.angle_internal_friction,
                              pressure,
                              partial_strain_rates[3],
                              p.max_yield_stress);
              }
            else
              {
                drpr_stress = creep_stress;
              }

            if ((std::fabs((diff_stress - creep_stress)/creep_stress) > 1e-6)
                || (std::fabs((disl_stress - creep_stress)/creep_stress) > 1e-6)
                || (std::fabs((prls_stress - creep_stress)/creep_stress) > 1e-6)
                || (std::fabs((drpr_stress - creep_stress)/creep_stress) > 1e-6))
              {
                error = true;
                std::cout << "   creep stress: " << creep_stress;
                std::cout << " diffusion stress: " << diff_stress;
                std::cout << " dislocation stress: " << disl_stress;
                std::cout << " peierls stress: " << prls_stress;
                std::cout << " drucker prager stress: " << drpr_stress << std::endl;
              }
          }
      }

    if (error)
      {
        std::cout << "   Error: The individual creep stresses differ by more than the required tolerance." << std::endl;
        std::cout << "Some parts of the test were not successful." << std::endl;
      }
    else
      {
        std::cout << "OK" << std::endl;
      }
  }

  template <>
  void f(const aspect::SimulatorAccess<2> &,
         aspect::Assemblers::Manager<2> &)
  {
    AssertThrow(false,dealii::ExcInternalError());
  }

  template <int dim>
  void signal_connector (aspect::SimulatorSignals<dim> &signals)
  {
    std::cout << "* Connecting signals" << std::endl;
    signals.set_assemblers.connect (std::bind(&f<dim>,
                                              std::placeholders::_1,
                                              std::placeholders::_2));
  }

  ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                    signal_connector<3>)
}
