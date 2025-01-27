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

      AssertThrow(material_file_names.size() == 1 || SimulatorAccess<dim>::get_end_time () == 0,
                  ExcMessage("The 'entropy model' material model can only handle one composition, "
                             "and can therefore only read one material lookup table."));



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
    void
    EntropyModel<dim>::evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                                MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      const unsigned int projected_density_index = this->introspection().compositional_index_for_name("density_field");
      //TODO : need to make it work for more than one field
      const std::vector<unsigned int> &entropy_indices = this->introspection().get_indices_for_fields_of_type(CompositionalFieldDescription::entropy);
      const unsigned int entropy_index = entropy_indices[0];
      const std::vector<unsigned int> &composition_indices = this->introspection().get_indices_for_fields_of_type(CompositionalFieldDescription::chemical_composition);

      AssertThrow(composition_indices.size() == material_file_names.size() - 1,
                  ExcMessage("The 'entropy model' material model assumes that there exists a background field in addition to the compositional fields, "
                             "and therefore it requires one more lookup table than there are chemical compositional fields."));

      EquationOfStateOutputs<dim> eos_outputs (material_file_names.size());
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
          const double entropy = in.composition[i][entropy_index];
          const double pressure = this->get_adiabatic_conditions().pressure(in.position[i]) / 1.e5;
          adjusted_inputs.temperature[i] = entropy_reader[0]->temperature(entropy,pressure);

          // Loop over all material files, and store the looked-up values for all compositions.
          for (unsigned int j=0; j<material_file_names.size(); ++j)
            {
              eos_outputs.densities[j] = entropy_reader[j]->density(entropy, pressure);
              eos_outputs.thermal_expansion_coefficients[j] = entropy_reader[j]->thermal_expansivity(entropy,pressure);
              eos_outputs.specific_heat_capacities[j] = entropy_reader[j]->specific_heat(entropy,pressure);

              const Tensor<1, 2> pressure_unit_vector({0.0, 1.0});
              eos_outputs.compressibilities[j] = ((entropy_reader[j]->density_gradient(entropy,pressure)) * pressure_unit_vector) / eos_outputs.densities[j];
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

          out.entropy_derivative_pressure[i]    = 0.;
          out.entropy_derivative_temperature[i] = 0.;
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c]            = 0.;

          // set up variable to interpolate prescribed field outputs onto compositional fields
          if (PrescribedFieldOutputs<dim> *prescribed_field_out = out.template get_additional_output<PrescribedFieldOutputs<dim>>())
            {
              prescribed_field_out->prescribed_field_outputs[i][projected_density_index] = out.densities[i];
            }

          // set up variable to interpolate prescribed field outputs onto temperature field
          if (PrescribedTemperatureOutputs<dim> *prescribed_temperature_out = out.template get_additional_output<PrescribedTemperatureOutputs<dim>>())
            {
              prescribed_temperature_out->prescribed_temperature_outputs[i] = adjusted_inputs.temperature[i];
            }

          // Calculate Viscosity
          if (in.requests_property(MaterialProperties::viscosity))
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

                  const double strain_rate_effective = std::fabs(second_invariant(deviator(in.strain_rate[i])));

                  if (std::sqrt(strain_rate_effective) >= std::numeric_limits<double>::min())
                    {
                      const double pressure =  this->get_adiabatic_conditions().pressure(in.position[i]);
                      const double eta_plastic = drucker_prager_plasticity.compute_viscosity(cohesion,
                                                                                             angle_of_internal_friction,
                                                                                             pressure,
                                                                                             std::sqrt(strain_rate_effective),
                                                                                             std::numeric_limits<double>::infinity());

                      effective_viscosity = 1.0 / ( ( 1.0 /  eta_plastic  ) + ( 1.0 / (vis_lateral * viscosity_profile) ) );

                      PlasticAdditionalOutputs<dim> *plastic_out = out.template get_additional_output<PlasticAdditionalOutputs<dim>>();
                      if (plastic_out != nullptr)
                        {
                          plastic_out->cohesions[i] = cohesion;
                          plastic_out->friction_angles[i] = angle_of_internal_friction;
                          plastic_out->yielding[i] = eta_plastic < (vis_lateral * viscosity_profile) ? 1 : 0;
                        }
                    }
                  out.viscosities[i] = std::max(std::min(effective_viscosity,max_eta),min_eta);
                }
            }

          // fill seismic velocities outputs if they exist
          if (SeismicAdditionalOutputs<dim> *seismic_out = out.template get_additional_output<SeismicAdditionalOutputs<dim>>())
            {

              std::vector<double> vp (material_file_names.size());
              std::vector<double> vs (material_file_names.size());
              for (unsigned int j=0; j<material_file_names.size(); ++j)
                {
                  vp[j] = entropy_reader[j]->seismic_vp(entropy,pressure);
                  vs[j] = entropy_reader[j]->seismic_vs(entropy,pressure);
                }
              seismic_out->vp[i] = MaterialUtilities::average_value (volume_fractions, vp, MaterialUtilities::arithmetic);
              seismic_out->vs[i] = MaterialUtilities::average_value (volume_fractions, vs, MaterialUtilities::arithmetic);
            }
        }

      // Evaluate thermal conductivity. This has to happen after
      // the evaluation of the equation of state and calculation of temperature.
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
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Entropy model");
        {
          data_directory              = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));
          material_file_names          = Utilities::split_string_list(prm.get ("Material file name"));
          lateral_viscosity_file_name  = prm.get ("Lateral viscosity file name");
          min_eta                     = prm.get_double ("Minimum viscosity");
          max_eta                     = prm.get_double ("Maximum viscosity");
          max_lateral_eta_variation    = prm.get_double ("Maximum lateral viscosity variation");

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

          angle_of_internal_friction = prm.get_double ("Angle of internal friction") * constants::degree_to_radians;
          cohesion = prm.get_double("Cohesion");

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
      if (out.template get_additional_output<SeismicAdditionalOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::SeismicAdditionalOutputs<dim>> (n_points));
        }

      if (out.template get_additional_output<PrescribedFieldOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::PrescribedFieldOutputs<dim>>
            (n_points, this->n_compositional_fields()));
        }

      if (out.template get_additional_output<PrescribedTemperatureOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::PrescribedTemperatureOutputs<dim>>
            (n_points));
        }

      if (out.template get_additional_output<PlasticAdditionalOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<PlasticAdditionalOutputs<dim>> (n_points));
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
