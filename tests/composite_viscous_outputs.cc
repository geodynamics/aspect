#include <aspect/simulator.h>
#include <aspect/material_model/rheology/composite_visco_plastic.h>

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
  auto n_phases = std::make_shared<std::vector<unsigned int>>(1); // 1 phase per composition
  const unsigned int composition = 0;

  // Next, we initialise instances of the composite rheology and
  // individual creep mechanisms.
  std::unique_ptr<Rheology::CompositeViscoPlastic<dim>> composite_creep;
  composite_creep = std_cxx14::make_unique<Rheology::CompositeViscoPlastic<dim>>();
  composite_creep->initialize_simulator (simulator_access.get_simulator());
  composite_creep->declare_parameters(prm);
  prm.set("Use plastic damper", "true");
  prm.set("Plastic damper viscosity", "1.e17");
  prm.set("Maximum yield stress", "5e8");
  composite_creep->parse_parameters(prm, n_phases);

  std::unique_ptr<Rheology::DiffusionCreep<dim>> diffusion_creep;
  diffusion_creep = std_cxx14::make_unique<Rheology::DiffusionCreep<dim>>();
  diffusion_creep->initialize_simulator (simulator_access.get_simulator());
  diffusion_creep->declare_parameters(prm);
  diffusion_creep->parse_parameters(prm, n_phases);

  std::unique_ptr<Rheology::DislocationCreep<dim>> dislocation_creep;
  dislocation_creep = std_cxx14::make_unique<Rheology::DislocationCreep<dim>>();
  dislocation_creep->initialize_simulator (simulator_access.get_simulator());
  dislocation_creep->declare_parameters(prm);
  dislocation_creep->parse_parameters(prm, n_phases);

  std::unique_ptr<Rheology::PeierlsCreep<dim>> peierls_creep;
  peierls_creep = std_cxx14::make_unique<Rheology::PeierlsCreep<dim>>();
  peierls_creep->initialize_simulator (simulator_access.get_simulator());
  peierls_creep->declare_parameters(prm);
  peierls_creep->parse_parameters(prm); // multiple phases not yet implemented for peierls

  std::unique_ptr<Rheology::DruckerPrager<dim>> drucker_prager;
  drucker_prager = std_cxx14::make_unique<Rheology::DruckerPrager<dim>>();
  drucker_prager->declare_parameters(prm);
  prm.set("Use plastic damper", "true");
  prm.set("Plastic damper viscosity", "1.e17");
  prm.set("Maximum yield stress", "5e8");
  Rheology::DruckerPragerParameters p = drucker_prager->parse_parameters(1, prm);

  // The creep components are arranged in series with each other.
  // This package of components is then arranged in parallel with
  // a strain rate limiter with a constant viscosity lim_visc.
  // The whole system is then arranged in series with a viscosity limiter with
  // viscosity max_visc.
  // lim_visc is equal to (min_visc*max_visc)/(max_visc - min_visc)
  double min_visc = aspect::Utilities::parse_map_to_double_array(prm.get("Minimum viscosities"),
                                                                 list_of_composition_names,
                                                                 true,
                                                                 "Minimum viscosities",
                                                                 true,
                                                                 n_phases)[composition];
  double max_visc = aspect::Utilities::parse_map_to_double_array(prm.get("Maximum viscosities"),
                                                                 list_of_composition_names,
                                                                 true,
                                                                 "Maximum viscosities",
                                                                 true,
                                                                 n_phases)[composition];
  double lim_visc = (min_visc*max_visc)/(max_visc - min_visc);

  // Assign values to the variables which will be passed to compute_viscosity
  // The test involves pure shear calculations at 1 GPa and variable temperature
  double temperature;
  const double pressure = 1.e9;
  SymmetricTensor<2,dim> strain_rate;
  strain_rate[0][0] = -1e-11;
  strain_rate[0][1] = 0.;
  strain_rate[1][1] = 1e-11;
  strain_rate[2][0] = 0.;
  strain_rate[2][1] = 0.;
  strain_rate[2][2] = 0.;

  std::cout << "temperature (K)   eta (Pas)   creep stress (Pa)   edot_ii (/s)   edot_ii fractions (diff, disl, prls, drpr, max)" << std::endl;

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
  std::vector<double> partial_strain_rates;

  for (unsigned int i=0; i <= 10; i++)
    {
      temperature = 1000. + i*100.;

      // Compute the viscosity
      viscosity = composite_creep->compute_viscosity(pressure, temperature, composition, strain_rate, partial_strain_rates);
      total_strain_rate = std::accumulate(partial_strain_rates.begin(), partial_strain_rates.end(), 0.);

      // The creep strain rate is calculated by subtracting the strain rate
      // of the max viscosity dashpot from the total strain rate
      // The creep stress is then calculated by subtracting the stress running
      // through the strain rate limiter from the total stress
      creep_strain_rate = total_strain_rate - partial_strain_rates[4];
      creep_stress = 2.*(viscosity*total_strain_rate - lim_visc*creep_strain_rate);

      // Print the output
      std::cout << temperature << " " << viscosity << " " << creep_stress << " " << total_strain_rate;
      for (unsigned int i=0; i < partial_strain_rates.size(); ++i)
        {
          std::cout << " " << partial_strain_rates[i]/total_strain_rate;
        }
      std::cout << std::endl;

      // The following lines test that each individual creep mechanism
      // experiences the same creep stress

      // Each creep mechanism should experience the same stress
      diff_stress = 2.*partial_strain_rates[0]*diffusion_creep->compute_viscosity(pressure, temperature, composition);
      disl_stress = 2.*partial_strain_rates[1]*dislocation_creep->compute_viscosity(partial_strain_rates[1], pressure, temperature, composition);
      prls_stress = 2.*partial_strain_rates[2]*peierls_creep->compute_viscosity(partial_strain_rates[2], pressure, temperature, composition);
      if (partial_strain_rates[3] > 0.)
        {
          drpr_stress = 2.*partial_strain_rates[3]*drucker_prager->compute_viscosity(p.cohesions[0],
                                                                                     p.angles_internal_friction[0],
                                                                                     pressure,
                                                                                     partial_strain_rates[3],
                                                                                     p.max_yield_stress,
                                                                                     p.use_plastic_damper,
                                                                                     p.damper_viscosity,
                                                                                     std::numeric_limits<double>::infinity());
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
  using namespace dealii;
  std::cout << "* Connecting signals" << std::endl;
  signals.set_assemblers.connect (std::bind(&f<dim>,
                                            std::placeholders::_1,
                                            std::placeholders::_2));
}

ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                  signal_connector<3>)
