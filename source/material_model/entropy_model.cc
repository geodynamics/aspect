/*
  Copyright (C) 2016 - 2024 by the authors of the ASPECT code.

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


#include <aspect/material_model/entropy_model.h>
#include <aspect/material_model/thermal_conductivity/constant.h>
#include <aspect/material_model/thermal_conductivity/tosi_stackhouse.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/utilities.h>

#include <deal.II/base/table.h>
#include <fstream>
#include <iostream>
#include <aspect/material_model/rheology/visco_plastic.h>
#include <aspect/material_model/steinberger.h>
#include <aspect/material_model/equation_of_state/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace
    {
      template <int dim>
      bool solver_scheme_is_supported(const Parameters<dim> &parameters)
      {
        // If we solve advection equations, we need to iterate them, because this material
        // models splits temperature diffusion from entropy advection.
        switch (parameters.nonlinear_solver)
          {
            case Parameters<dim>::NonlinearSolver::Kind::iterated_Advection_and_Stokes:
            case Parameters<dim>::NonlinearSolver::Kind::iterated_Advection_and_defect_correction_Stokes:
            case Parameters<dim>::NonlinearSolver::Kind::iterated_Advection_and_Newton_Stokes:
            case Parameters<dim>::NonlinearSolver::Kind::no_Advection_no_Stokes:
            case Parameters<dim>::NonlinearSolver::Kind::no_Advection_iterated_defect_correction_Stokes:
            case Parameters<dim>::NonlinearSolver::Kind::no_Advection_iterated_Stokes:
            case Parameters<dim>::NonlinearSolver::Kind::no_Advection_single_Stokes:
            case Parameters<dim>::NonlinearSolver::Kind::first_timestep_only_single_Stokes:
              return true;

            case Parameters<dim>::NonlinearSolver::Kind::single_Advection_single_Stokes:
            case Parameters<dim>::NonlinearSolver::Kind::single_Advection_iterated_Stokes:
            case Parameters<dim>::NonlinearSolver::Kind::single_Advection_iterated_defect_correction_Stokes:
            case Parameters<dim>::NonlinearSolver::Kind::single_Advection_iterated_Newton_Stokes:
            case Parameters<dim>::NonlinearSolver::Kind::single_Advection_no_Stokes:
              return false;
          }
        Assert(false, ExcNotImplemented());
        return false;
      }
    }



    template <int dim>
    void
    EntropyModel<dim>::initialize()
    {
      CitationInfo::add("entropy");

      AssertThrow (this->get_parameters().formulation_mass_conservation ==
                   Parameters<dim>::Formulation::MassConservation::projected_density_field,
                   ExcMessage("The 'entropy model' material model was only tested with the "
                              "'projected density field' approximation "
                              "for the mass conservation equation, which is not selected."));

      AssertThrow (this->introspection().composition_type_exists(CompositionalFieldDescription::Type::entropy),
                   ExcMessage("The 'entropy model' material model requires the existence of a compositional field "
                              "named 'entropy'. This field does not exist."));

      AssertThrow(solver_scheme_is_supported(this->get_parameters()) == true,
                  ExcMessage("The 'entropy model' material model requires the use of a solver scheme that "
                             "iterates over the advection equations but a non iterating solver scheme was selected. "
                             "Please check the consistency of your solver scheme."));

      for (unsigned int i = 0; i < material_file_names.size(); ++i)
        {
          entropy_reader.push_back(std::make_unique<MaterialUtilities::Lookup::EntropyReader>());
          entropy_reader[i]->initialize(this->get_mpi_communicator(), data_directory, material_file_names[i]);
        }

      lateral_viscosity_prefactor_lookup = std::make_unique<internal::LateralViscosityLookup>(data_directory+lateral_viscosity_file_name,
                                           this->get_mpi_communicator());
    }



    template <int dim>
    bool
    EntropyModel<dim>::
    is_compressible () const
    {
      return true;
    }



    template <int dim>
    double
    EntropyModel<dim>::
    equilibrate_temperature (const std::vector<double> &temperature,
                             const std::vector<double> &mass_fractions,
                             const std::vector<double> &entropy,
                             const std::vector<double> &specific_heat,
                             const double pressure,
                             std::vector<double> &component_equilibrated_S
                            ) const
    {
      AssertThrow(material_file_names.size() == temperature.size() &&
                  temperature.size() == mass_fractions.size() &&
                  temperature.size() == entropy.size() &&
                  temperature.size() == specific_heat.size(),
                  ExcMessage("The temperature, chemical composition, entropy and specific heat capacity vectors"
                             " must all have the same size as the number of look-up tables."));

      std::vector<double> component_entropies = entropy;
      std::vector<double> component_temperatures = temperature;
      std::vector<double> component_heat_capacities = specific_heat;

      bool equilibration = false;
      unsigned int iteration = 0;
      double ln_equilibrated_T = 0;

      // We do the following thermodynamic iteration to find the equilibrated temperature of our components
      // while conserving the total entropy.
      while (equilibration == false)
        {
          if (iteration >= multicomponent_max_iteration)
            {
              std::ostringstream error_message;
              error_message << "The components are not equilibrated after the maximum number of iterations: " << std::to_string(multicomponent_max_iteration)
                            << ". Equilibration failed at pressure = " << std::to_string(pressure) << ". The last equilibrated T = " << std::to_string(std::exp(ln_equilibrated_T))
                            << ". Individual component temperatures:";

              for (unsigned int k=0; k<component_temperatures.size(); ++k)
                {
                  error_message << " Component temperature [" << std::to_string(k) << "] = " << std::to_string(component_temperatures[k]) << '.';
                }

              Assert(false,
                     ExcMessage(error_message.str()));

              std::cerr << error_message.str() << std::endl;
              MPI_Abort(MPI_COMM_WORLD, 1);
            }

          double T_numerator = 0;
          double T_denominator = 0;

          // Step 1: Guess the equilibrated temperature based on the components' entropies
          for (unsigned int i = 0; i < material_file_names.size(); ++i)
            {
              T_numerator += mass_fractions[i] * component_heat_capacities[i] * std::log(component_temperatures[i]);
              T_denominator += mass_fractions[i] * component_heat_capacities[i];
            }

          ln_equilibrated_T = T_numerator/T_denominator;

          // Step 2: Calculate the new entropies for each component based on the guessed equilibrated temperature.
          // The components' entropies are updated with the iteration
          for (unsigned int i = 0; i < material_file_names.size(); ++i)
            {
              component_entropies[i] = component_entropies[i] + component_heat_capacities[i] * (ln_equilibrated_T - std::log (component_temperatures[i]));

              // Step 3 Look up the new temperature and specific heat capacity for each component based on the new entropies
              // The components' entropies are updated with the iteration.
              component_temperatures[i] = entropy_reader[i]->temperature(component_entropies[i], pressure);
              component_heat_capacities[i] = entropy_reader[i]->specific_heat(component_entropies[i], pressure);
            }

          equilibration = true;

          // Step 4: Are the new temperatures equal to the guessed equilibrated temperature?
          // If not, we need to iterate again
          // If yes, we return the equilibrated temperature
          for (unsigned int i = 0; i < material_file_names.size(); ++i)
            {
              // If a component has mass fraction of 0, we do not check its temperature.
              // It is because its entropy and other properties does not affect the equilibrated temperature calculation.
              // Therefore, its temperature is meaningless and does not need to be equilibrated.
              // For example, if there are 2 components including background, and one of them has 0 mass fraction,
              // then the equilibrated temperature is always the same as the 100% mass fraction component's temperature.
              if (mass_fractions[i] > 0. && std::abs (component_temperatures[i] - std::exp(ln_equilibrated_T)) >= multicomponent_tolerance)
                {
                  equilibration = false;
                  break;
                }
            }

          ++iteration;
        }

      component_equilibrated_S = component_entropies;
      return std::exp(ln_equilibrated_T);
    }



    template <int dim>
    void
    EntropyModel<dim>::evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                                MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      const unsigned int projected_density_index = this->introspection().compositional_index_for_name("density_field");

      const std::vector<unsigned int> &entropy_indices = this->introspection().get_indices_for_fields_of_type(CompositionalFieldDescription::entropy);
      const std::vector<unsigned int> &composition_indices = this->introspection().get_indices_for_fields_of_type(CompositionalFieldDescription::chemical_composition);

      AssertThrow(composition_indices.size() == material_file_names.size() - 1,
                  ExcMessage("The 'entropy model' material model assumes that there exists a background field in addition to the compositional fields, "
                             "and therefore it requires one more lookup table than there are chemical compositional fields."));

      EquationOfStateOutputs<dim> eos_outputs (material_file_names.size());
      const std::shared_ptr<ReactionRateOutputs<dim>> reaction_rate_out = out.template get_additional_output_object<ReactionRateOutputs<dim>>();
      std::vector<double> volume_fractions (material_file_names.size());
      std::vector<double> mass_fractions (material_file_names.size());

      // We need to make a copy of the material model inputs because we want to replace the
      // temperature with the temperature from the lookup table.
      MaterialModel::MaterialModelInputs<dim> adjusted_inputs(in);

      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          // Use the adiabatic pressure instead of the real one,
          // to stabilize against pressure oscillations in phase transitions.
          // This is a requirement of the projected density approximation for
          // the Stokes equation and not related to the entropy formulation.
          // Also convert pressure from Pa to bar, bar is used in the table.

          std::vector<double> component_entropy (material_file_names.size());
          std::vector<double> component_temperature_lookup (material_file_names.size());
          const double pressure = this->get_adiabatic_conditions().pressure(in.position[i]) / 1.e5;

          // Loop over all material files, and store the looked-up values for all components.
          for (unsigned int j=0; j<material_file_names.size(); ++j)
            {
              component_entropy[j] = in.composition[i][entropy_indices[j]];
              component_temperature_lookup[j] = entropy_reader[j]->temperature(component_entropy[j], pressure);

              eos_outputs.densities[j] = entropy_reader[j]->density(component_entropy[j], pressure);
              eos_outputs.thermal_expansion_coefficients[j] = entropy_reader[j]->thermal_expansivity(component_entropy[j],pressure);
              eos_outputs.specific_heat_capacities[j] = entropy_reader[j]->specific_heat(component_entropy[j],pressure);

              const Tensor<1, 2> pressure_unit_vector({0.0, 1.0});
              eos_outputs.compressibilities[j] = ((entropy_reader[j]->density_gradient(component_entropy[j],pressure)) * pressure_unit_vector) / eos_outputs.densities[j];
            }

          // Calculate volume fractions from mass fractions
          // If there is only one lookup table, set the mass and volume fractions to 1
          if (material_file_names.size() == 1)
            mass_fractions [0] = 1.0;

          else
            {
              // We only want to compute mass/volume fractions for fields that are chemical compositions.
              mass_fractions = MaterialUtilities::compute_only_composition_fractions(in.composition[i], this->introspection().chemical_composition_field_indices());
            }

          volume_fractions = MaterialUtilities::compute_volumes_from_masses(mass_fractions,
                                                                            eos_outputs.densities,
                                                                            true);

          out.densities[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.densities, MaterialUtilities::arithmetic);
          out.thermal_expansion_coefficients[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.thermal_expansion_coefficients, MaterialUtilities::arithmetic);

          out.specific_heat[i] = MaterialUtilities::average_value (mass_fractions, eos_outputs.specific_heat_capacities, MaterialUtilities::arithmetic);
          out.compressibilities[i] = MaterialUtilities::average_value (mass_fractions, eos_outputs.compressibilities, MaterialUtilities::arithmetic);

          // The component_equilibrated_S will be updated while we are looking for the equilibrated temperature.
          // It is used to calculate the reaction terms later.
          std::vector<double> component_equilibrated_S (material_file_names.size());

          const double equilibrated_T = equilibrate_temperature (component_temperature_lookup,
                                                                 mass_fractions, component_entropy,
                                                                 eos_outputs.specific_heat_capacities,
                                                                 pressure,
                                                                 component_equilibrated_S);

          adjusted_inputs.temperature[i] = equilibrated_T;

          out.entropy_derivative_pressure[i]    = 0.;
          out.entropy_derivative_temperature[i] = 0.;

          // Calculate the reaction terms
          if (material_file_names.size() == 1)
            {
              for (unsigned int c=0; c<in.composition[i].size(); ++c)
                {
                  out.reaction_terms[i][c] = 0.;
                }
            }

          else
            {
              // Calculate the reaction rates for the operator splitting
              for (unsigned int c = 0; c < in.composition[i].size(); ++c)
                {
                  if (this->get_parameters().use_operator_splitting && reaction_rate_out != nullptr)
                    {
                      reaction_rate_out->reaction_rates[i][c] = 0.0;

                      // Figure out if compositional field c is an entropy field and the how manyth entropy field it is
                      bool c_is_entropy_field = false;
                      unsigned int c_is_nth_entropy_field = 0;

                      unsigned int nth_entropy_index = 0;
                      for (unsigned int entropy_index : entropy_indices)
                        {
                          if (c == entropy_index)
                            {
                              c_is_entropy_field = true;
                              c_is_nth_entropy_field = nth_entropy_index;
                            }
                          ++nth_entropy_index;
                        }

                      const unsigned int timestep_number = this->simulator_is_past_initialization()
                                                           ?
                                                           this->get_timestep_number()
                                                           :
                                                           0;

                      if (c_is_entropy_field == true && timestep_number > 0)
                        reaction_rate_out->reaction_rates[i][c] = (component_equilibrated_S[c_is_nth_entropy_field] - in.composition[i][entropy_indices[c_is_nth_entropy_field]]) / this->get_timestep();
                    }

                  out.reaction_terms[i][c] = 0.0;
                }
            }

          // set up variable to interpolate prescribed field outputs onto compositional fields
          if (const std::shared_ptr<PrescribedFieldOutputs<dim>> prescribed_field_out
              = out.template get_additional_output_object<PrescribedFieldOutputs<dim>>())
            {
              prescribed_field_out->prescribed_field_outputs[i][projected_density_index] = out.densities[i];
            }

          // set up variable to interpolate prescribed field outputs onto temperature field
          if (const std::shared_ptr<PrescribedTemperatureOutputs<dim>> prescribed_temperature_out
              = out.template get_additional_output_object<PrescribedTemperatureOutputs<dim>>())
            {
              prescribed_temperature_out->prescribed_temperature_outputs[i] = adjusted_inputs.temperature[i];
            }

          // Calculate Viscosity
          if (in.requests_property(MaterialProperties::viscosity) || in.requests_property(MaterialProperties::additional_outputs))
            {
              // read in the viscosity profile
              const double depth = this->get_geometry_model().depth(in.position[i]);
              const double viscosity_profile = depth_dependent_rheology->compute_viscosity(depth);

              // lateral viscosity variations
              const double reference_temperature = this->get_adiabatic_conditions().is_initialized()
                                                   ?
                                                   this->get_adiabatic_conditions().temperature(in.position[i])
                                                   :
                                                   this->get_parameters().adiabatic_surface_temperature;

              const double delta_temperature = adjusted_inputs.temperature[i] - reference_temperature;

              // Steinberger & Calderwood viscosity
              if (adjusted_inputs.temperature[i]*reference_temperature == 0)
                out.viscosities[i] = max_eta;
              else
                {
                  double vis_lateral = std::exp(-lateral_viscosity_prefactor_lookup->lateral_viscosity(depth)*delta_temperature/(adjusted_inputs.temperature[i]*reference_temperature));
                  // lateral vis variation
                  vis_lateral = std::max(std::min((vis_lateral),max_lateral_eta_variation),1/max_lateral_eta_variation);

                  if (std::isnan(vis_lateral))
                    vis_lateral = 1.0;

                  double effective_viscosity = vis_lateral * viscosity_profile;

                  const double pressure = this->get_adiabatic_conditions().pressure(in.position[i]);

                  MaterialModel::Rheology::DruckerPragerParameters drucker_prager_parameters;
                  drucker_prager_parameters.cohesion = cohesion;
                  drucker_prager_parameters.angle_internal_friction = angle_of_internal_friction;
                  drucker_prager_parameters.max_yield_stress = std::numeric_limits<double>::infinity();

                  const std::shared_ptr<PlasticAdditionalOutputs<dim>>
                  plastic_out = out.template get_additional_output_object<PlasticAdditionalOutputs<dim>>();

                  if (plastic_out != nullptr && in.requests_property(MaterialProperties::additional_outputs))
                    {
                      plastic_out->cohesions[i] = cohesion;
                      plastic_out->friction_angles[i] = angle_of_internal_friction;
                      plastic_out->yielding[i] = 0;
                      plastic_out->yield_stresses[i] = drucker_prager_plasticity.compute_yield_stress(pressure,
                                                                                                      drucker_prager_parameters);
                    }

                  const double strain_rate_effective = std::fabs(second_invariant(deviator(in.strain_rate[i])));

                  if (std::sqrt(strain_rate_effective) >= std::numeric_limits<double>::min())
                    {
                      const double eta_plastic = drucker_prager_plasticity.compute_viscosity(pressure,
                                                                                             std::sqrt(strain_rate_effective),
                                                                                             drucker_prager_parameters);

                      effective_viscosity = 1.0 / ( ( 1.0 /  eta_plastic  ) + ( 1.0 / (vis_lateral * viscosity_profile) ) );

                      if (plastic_out != nullptr && in.requests_property(MaterialProperties::additional_outputs))
                        plastic_out->yielding[i] = eta_plastic < (vis_lateral * viscosity_profile) ? 1 : 0;
                    }

                  out.viscosities[i] = std::max(std::min(effective_viscosity,max_eta),min_eta);
                }
            }

          // fill seismic velocities outputs if they exist
          if (const std::shared_ptr<SeismicAdditionalOutputs<dim>> seismic_out
              = out.template get_additional_output_object<SeismicAdditionalOutputs<dim>>())
            if (in.requests_property(MaterialProperties::additional_outputs))
              {
                std::vector<double> vp (material_file_names.size());
                std::vector<double> vs (material_file_names.size());
                for (unsigned int j=0; j<material_file_names.size(); ++j)
                  {
                    vp[j] = entropy_reader[j]->seismic_vp(component_entropy[j],pressure);
                    vs[j] = entropy_reader[j]->seismic_vs(component_entropy[j],pressure);
                  }

                seismic_out->vp[i] = MaterialUtilities::average_value (volume_fractions, vp, MaterialUtilities::arithmetic);
                seismic_out->vs[i] = MaterialUtilities::average_value (volume_fractions, vs, MaterialUtilities::arithmetic);
              }
        }

      // Evaluate thermal conductivity. This has to happen after
      // the evaluation of the equation of state and calculation of temperature.
      // Thermal conductivity can be pressure temperature dependent
      thermal_conductivity->evaluate(adjusted_inputs, out);
    }



    template <int dim>
    void
    EntropyModel<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Entropy model");
        {
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/material-model/entropy-table/opxtable/",
                             Patterns::DirectoryName (),
                             "The path to the model data. The path may also include the special "
                             "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the `data/' subdirectory of ASPECT.");
          prm.declare_entry ("Material file name", "material_table.txt",
                             Patterns::List (Patterns::Anything()),
                             "The file name of the material data. The first material data file is intended for the background composition. ");
          prm.declare_entry ("Reference viscosity", "1e22",
                             Patterns::Double(0),
                             "The viscosity that is used in this model. "
                             "\n\n"
                             "Units: \\si{\\pascal\\second}");
          prm.declare_entry ("Lateral viscosity file name", "temp-viscosity-prefactor.txt",
                             Patterns::Anything (),
                             "The file name of the lateral viscosity prefactor.");
          prm.declare_entry ("Minimum viscosity", "1e19",
                             Patterns::Double (0.),
                             "The minimum viscosity that is allowed in the viscosity "
                             "calculation. Smaller values will be cut off.");
          prm.declare_entry ("Maximum viscosity", "1e23",
                             Patterns::Double (0.),
                             "The maximum viscosity that is allowed in the viscosity "
                             "calculation. Larger values will be cut off.");
          prm.declare_entry ("Maximum lateral viscosity variation", "1e2",
                             Patterns::Double (0.),
                             "The relative cutoff value for lateral viscosity variations "
                             "caused by temperature deviations. The viscosity may vary "
                             "laterally by this factor squared.");
          prm.declare_entry ("Angle of internal friction", "0.",
                             Patterns::Double (0.),
                             "The value of the angle of internal friction, $\\phi$."
                             "For a value of zero, in 2D the von Mises criterion is retrieved. "
                             "Angles higher than 30 degrees are harder to solve numerically."
                             "Units: degrees.");
          prm.declare_entry ("Cohesion", "1e20",
                             Patterns::Double (0.),
                             "The value of the cohesion, $C$. The extremely large default"
                             "cohesion value (1e20 Pa) prevents the viscous stress from "
                             "exceeding the yield stress. Units: \\si{\\pascal}.");

          // Multicomponent equilibration parameters
          prm.declare_entry ("Maximum iteration for multicomponent equilibration", "50",
                             Patterns::Double (0.),
                             "The maximum allowed number of iterations for the multicomponent equlibration "
                             "to reach the tolerance value. If the maximum iteration is reached but the"
                             "temperature has not been equilibrated, the model run will abort.");
          prm.declare_entry ("Multicomponent equilibration tolerance", "1e-7",
                             Patterns::Double (0.),
                             "This is the maximum temperature difference between the different compositions "
                             "when they are considered in equilibrium."
                            );

          // Thermal conductivity parameters
          ThermalConductivity::Constant<dim>::declare_parameters(prm);
          prm.declare_entry ("Thermal conductivity formulation", "constant",
                             Patterns::Selection("constant|p-T-dependent"),
                             "Which law should be used to compute the thermal conductivity. "
                             "The 'constant' law uses a constant value for the thermal "
                             "conductivity. The 'p-T-dependent' formulation uses equations "
                             "from Stackhouse et al. (2015): First-principles calculations "
                             "of the lattice thermal conductivity of the lower mantle "
                             "(https://doi.org/10.1016/j.epsl.2015.06.050), and Tosi et al. "
                             "(2013): Mantle dynamics with pressure- and temperature-dependent "
                             "thermal expansivity and conductivity "
                             "(https://doi.org/10.1016/j.pepi.2013.02.004) to compute the "
                             "thermal conductivity in dependence of temperature and pressure. "
                             "The thermal conductivity parameter sets can be chosen in such a "
                             "way that either the Stackhouse or the Tosi relations are used. "
                             "The conductivity description can consist of several layers with "
                             "different sets of parameters. Note that the Stackhouse "
                             "parametrization is only valid for the lower mantle (bridgmanite).");
          ThermalConductivity::TosiStackhouse<dim>::declare_parameters(prm);

          prm.leave_subsection();
        }

        // Depth-dependent parameters from the rheology plugin
        Rheology::AsciiDepthProfile<dim>::declare_parameters(prm,
                                                             "Depth dependent viscosity");

        prm.leave_subsection();
      }
    }



    template <int dim>
    void
    EntropyModel<dim>::parse_parameters (ParameterHandler &prm)
    {
      if (this->introspection().composition_type_exists(CompositionalFieldDescription::Type::chemical_composition))
        {
          AssertThrow (this->get_parameters().use_operator_splitting == true,
                       ExcMessage("The 'entropy model' material model requires the use of operator splitting for multiple chemical composition. "
                                  "Please set the 'Use operator splitting' parameter to 1 in the material model parameters."));
        }

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Entropy model");
        {
          data_directory               = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));
          material_file_names          = Utilities::split_string_list(prm.get ("Material file name"));

          // Viscosity parameters
          lateral_viscosity_file_name  = prm.get ("Lateral viscosity file name");
          min_eta                      = prm.get_double ("Minimum viscosity");
          max_eta                      = prm.get_double ("Maximum viscosity");
          max_lateral_eta_variation    = prm.get_double ("Maximum lateral viscosity variation");

          // Placiticity parameters
          angle_of_internal_friction   = prm.get_double ("Angle of internal friction") * constants::degree_to_radians;
          cohesion                     = prm.get_double("Cohesion");

          // Multicomponent equilibration parameters
          multicomponent_max_iteration = prm.get_double("Maximum iteration for multicomponent equilibration");
          multicomponent_tolerance = prm.get_double("Multicomponent equilibration tolerance");

          // Thermal conductivity parameters
          if (prm.get ("Thermal conductivity formulation") == "constant")
            thermal_conductivity = std::make_unique<ThermalConductivity::Constant<dim>>();
          else if (prm.get ("Thermal conductivity formulation") == "p-T-dependent")
            {
              thermal_conductivity = std::make_unique<ThermalConductivity::TosiStackhouse<dim>>();
              if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(thermal_conductivity.get()))
                sim->initialize_simulator (this->get_simulator());
            }
          else
            AssertThrow(false, ExcMessage("Not a valid thermal conductivity formulation"));

          thermal_conductivity->parse_parameters(prm);

          prm.leave_subsection();
        }

        depth_dependent_rheology = std::make_unique<Rheology::AsciiDepthProfile<dim>>();
        depth_dependent_rheology->initialize_simulator (this->get_simulator());
        depth_dependent_rheology->parse_parameters(prm, "Depth dependent viscosity");
        depth_dependent_rheology->initialize();

        prm.leave_subsection();

        // Declare dependencies on solution variables
        this->model_dependence.viscosity = NonlinearDependence::temperature;
        this->model_dependence.density = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
        this->model_dependence.compressibility = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
        this->model_dependence.specific_heat = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
        this->model_dependence.thermal_conductivity = NonlinearDependence::none;
      }
    }



    template <int dim>
    void
    EntropyModel<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template has_additional_output_object<SeismicAdditionalOutputs<dim>>() == false)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::SeismicAdditionalOutputs<dim>> (n_points));
        }

      if (out.template has_additional_output_object<PrescribedFieldOutputs<dim>>() == false)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::PrescribedFieldOutputs<dim>>
            (n_points, this->n_compositional_fields()));
        }

      if (out.template has_additional_output_object<PrescribedTemperatureOutputs<dim>>() == false)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::PrescribedTemperatureOutputs<dim>>
            (n_points));
        }

      if (out.template has_additional_output_object<PlasticAdditionalOutputs<dim>>() == false)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<PlasticAdditionalOutputs<dim>> (n_points));
        }

      if (this->get_parameters().use_operator_splitting
          && out.template get_additional_output_object<ReactionRateOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::ReactionRateOutputs<dim>> (n_points, this->n_compositional_fields()));
        }
    }
  }
}



// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(EntropyModel,
                                   "entropy model",
                                   "A material model that is designed to use pressure and entropy (rather "
                                   "than pressure and temperature) as independent variables. "
                                   "It requires a thermodynamic data table that contains "
                                   "all relevant properties in a specific format as illustrated in "
                                   "the data/material-model/entropy-table/opxtable example folder. "
                                   "The material model requires the use of the projected density "
                                   "approximation for compressibility, and the existence of a "
                                   "compositional field called 'entropy'.")
  }
}
