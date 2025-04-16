/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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

#include <aspect/material_model/visco_plastic_simple.h>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    ViscoPlasticSimple<dim>::
    calculate_average_viscosity(const std::vector<std::vector<double>> &volume_fractions,
                                const double strain_rate,
                                const double pressure,
                                const double temperature,
                                const double plastic_strain,
                                const bool fill_reaction_terms_for_plastic_strain,
                                MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      const unsigned int n_evaluation_points = volume_fractions.size();
      const unsigned int n_chemical_fields = volume_fractions[0].size();

      std::vector<double> composition_viscosities(n_chemical_fields, numbers::signaling_nan<double>());
      std::vector<double> composition_plastic_strain_rates(n_chemical_fields, numbers::signaling_nan<double>());

      for (unsigned int j = 0; j < n_chemical_fields; ++j)
        {
          // Do not calculate the viscosity for the current compositional field if
          // its volume fraction is zero
          double max_fraction = 0.;
          for (unsigned int i = 0; i < n_evaluation_points; ++i)
            max_fraction = std::max(max_fraction, volume_fractions[i][j]);
          if (max_fraction < 1.e-8)
            {
              composition_viscosities[j] = max_viscosity;
              continue;
            }

          // Step 1: calculate the non-yielding viscosity
          double non_yielding_stress    = numbers::signaling_nan<double>();
          double non_yielding_viscosity = numbers::signaling_nan<double>();

          const Rheology::DislocationCreepParameters dislocation_creep_parameters = dislocation_creep.compute_creep_parameters(j);
          const Rheology::DiffusionCreepParameters diffusion_creep_parameters = diffusion_creep.compute_creep_parameters(j);

          // Check if the non-yielding viscosity can be calculated straightforwardly.
          // It occurs under two situations:
          // (a) The viscoelastic rheology is governed by a single constitutive relation;
          // (b) The viscoelastic rheology is linear.
          const bool single_constitutive_relation = (this->get_parameters().enable_elasticity == false && viscous_flow_law != composite);
          const bool rheology_is_nonlinear = (viscous_flow_law != dislocation && diffusion_creep_parameters.stress_exponent != 1.) ||
                                             (viscous_flow_law != diffusion && dislocation_creep_parameters.stress_exponent != 1.);

          if (single_constitutive_relation == true || rheology_is_nonlinear == false)
            {
              switch (viscous_flow_law)
                {
                  case diffusion:
                  {
                    non_yielding_viscosity = diffusion_creep.compute_viscosity(pressure, temperature, j);
                    break;
                  }
                  case dislocation:
                  {
                    non_yielding_viscosity = dislocation_creep.compute_viscosity(strain_rate, pressure, temperature, j);
                    break;
                  }

                  case composite:
                  {
                    const double viscosity_diffusion   = diffusion_creep.compute_viscosity(pressure, temperature, j);
                    const double viscosity_dislocation = dislocation_creep.compute_viscosity(strain_rate, pressure, temperature, j);
                    non_yielding_viscosity = (viscosity_diffusion * viscosity_dislocation) /
                                             (viscosity_diffusion + viscosity_dislocation);
                    break;
                  }

                  default:
                  {
                    AssertThrow(false, ExcNotImplemented());
                  }
                }

              if (this->get_parameters().enable_elasticity)
                non_yielding_viscosity = elastic_rheology.calculate_viscoelastic_viscosity(non_yielding_viscosity,
                                                                                           elastic_rheology.get_elastic_shear_moduli()[j]);

              non_yielding_stress = 2. * non_yielding_viscosity * strain_rate;
            }
          else
            {
              // If the viscoelastic rheology is governed by multiple constitutive relations and
              // at least one of them is nonlinear, then we solve the second invariant of deviatoric
              // stress with Newton's method.
              unsigned int iteration  = 0;
              double initial_residual = numbers::signaling_nan<double>();
              double residual         = numbers::signaling_nan<double>();

              // Start with the assumption that the strain rate is accomodated by one viscous flow law.
              double viscosity_initial_guess = 
                (viscous_flow_law != dislocation ?
                 diffusion_creep.compute_viscosity(pressure, temperature, j) :
                 dislocation_creep.compute_viscosity(strain_rate, pressure, temperature, j));

              double log_tau = std::log(2. * viscosity_initial_guess * strain_rate);

              do
                {
                  // r  = log(edot_disl + edot_diff + edot_e) - log_edot_star
                  // r' = (edot_disl * dedot_disl/dtau + edot_diff * dedot_diff/dtau + edot_e * dedot_e/dtau)
                  //      / (edot_disl + edot_diff + edot_e)
                  const double log_edot_star = std::log(strain_rate);
                  
                  double numerator   = 0.;
                  double denominator = 0.;

                  if (viscous_flow_law != dislocation)
                    {
                      const std::pair<double, double> log_diff_edot_and_deriv =
                        diffusion_creep.compute_log_strain_rate_and_derivative(log_tau, pressure, temperature, diffusion_creep_parameters);

                      const double strain_rate_diffusion = std::exp(log_diff_edot_and_deriv.first);

                      numerator   += strain_rate_diffusion * log_diff_edot_and_deriv.second;
                      denominator += strain_rate_diffusion;
                    }

                  if (viscous_flow_law != diffusion)
                    {
                      const std::pair<double, double> log_disl_edot_and_deriv = 
                        dislocation_creep.compute_log_strain_rate_and_derivative(log_tau, pressure, temperature, dislocation_creep_parameters);

                      const double strain_rate_dislocation = std::exp(log_disl_edot_and_deriv.first);

                      numerator   += strain_rate_dislocation * log_disl_edot_and_deriv.second;
                      denominator += strain_rate_dislocation;
                    }

                  if (this->get_parameters().enable_elasticity)
                    {
                      const double strain_rate_elastic = std::exp(log_tau -
                        std::log(2. * elastic_rheology.calculate_elastic_viscosity(elastic_rheology.get_elastic_shear_moduli()[j])));

                      numerator   += strain_rate_elastic;
                      denominator += strain_rate_elastic;
                    }

                  AssertThrow(denominator > std::numeric_limits<double>::min(),
                              ExcMessage("The denominator in the expression of dr/dtau is smaller than or equal to zero."));

                  residual = std::log(denominator) - log_edot_star;
                  if (iteration == 0)
                    {
                      initial_residual = residual;
                      if (std::abs(initial_residual) < std::numeric_limits<double>::min())
                        break;
                    }

                  if (residual / initial_residual < newton_iteration_tolerance)
                    break;

                  const double dr_dlogtau = numerator / denominator;
                  AssertThrow(std::fabs(dr_dlogtau) > std::numeric_limits<double>::min(),
                              ExcMessage("The derivative of Newton residual is equal to zero."));

                  log_tau -= residual / dr_dlogtau;
                }
              while (iteration < max_newton_iteration);

              // Make sure that the iteration converges.
              AssertThrow(residual / initial_residual < newton_iteration_tolerance,
                          ExcMessage("The Newton iterations failed to converge!"));
              
              non_yielding_stress = std::exp(log_tau);
              non_yielding_viscosity = non_yielding_stress / (2. * strain_rate);
            }

          // Step 2: plastic yielding

          // Step 2a: Calculate the plastic strain weakening factors
          double weakening_cohesion = 1.;
          double weakening_friction = 1.;
          if (enable_plastic_strain_weakening)
            {
              const double cut_off_strain = std::max(std::min(plastic_strain,
                                                              end_plastic_strain_weakening_intervals[j]),
                                                     start_plastic_strain_weakening_intervals[j]);

              const double strain_fraction = (cut_off_strain - start_plastic_strain_weakening_intervals[j]) /
                                             (start_plastic_strain_weakening_intervals[j] - end_plastic_strain_weakening_intervals[j]);

              weakening_cohesion = 1. + (1. - cohesion_strain_weakening_factors[j]) * strain_fraction;
              weakening_friction = 1. + (1. - friction_strain_weakening_factors[j]) * strain_fraction;
            }

          // Step 2b: calculate the current friction and cohesion parameters
          const Rheology::DruckerPragerParameters drucker_prager_parameters = drucker_prager_plasticity.compute_drucker_prager_parameters(j);

          double current_cohesion = drucker_prager_parameters.cohesion * weakening_cohesion;
          double current_friction = drucker_prager_parameters.angle_internal_friction * weakening_friction;
          current_friction = friction_models.compute_friction_angle(strain_rate, j, current_friction, 
                                                                    /*dummy position*/Point<dim>());

          // Step 2c: calculate the Drucker-Prager yield stress
          double pressure_for_plasticity = pressure;
          if (allow_negative_pressure_in_plasticity == false)
            pressure_for_plasticity = std::max(pressure_for_plasticity, 0.);

          const double yield_stress = drucker_prager_plasticity.compute_yield_stress(current_cohesion,
                                                                                     current_friction,
                                                                                     pressure_for_plasticity,
                                                                                     drucker_prager_parameters.max_yield_stress);

          double effective_viscosity = non_yielding_viscosity;
          double plastic_strain_rate = 0.;
          if (non_yielding_stress > yield_stress)
            {
              effective_viscosity = drucker_prager_plasticity.compute_viscosity(current_cohesion,
                                                                                current_friction,
                                                                                pressure_for_plasticity,
                                                                                strain_rate,
                                                                                drucker_prager_parameters.max_yield_stress,
                                                                                non_yielding_viscosity);

              plastic_strain_rate = strain_rate - yield_stress / non_yielding_viscosity;
            }

          composition_viscosities[j] = std::min(std::max(effective_viscosity, min_viscosity), max_viscosity);
          composition_plastic_strain_rates[j] = plastic_strain_rate;
        }

      // First calculate the average viscosity at each evaluation point, then perform cellwise averaging
      std::vector<double> point_viscosities(n_evaluation_points);
      for (unsigned int i = 0; i < n_evaluation_points; ++i)
         point_viscosities[i] = MaterialUtilities::average_value(volume_fractions[i], composition_viscosities, pointwise_viscosity_averaging);

      const typename DoFHandler<dim>::active_cell_iterator current_cell(&this->get_triangulation(),
                                                                        fe_values->get_cell()->level(),
                                                                        fe_values->get_cell()->index(),
                                                                        &this->get_dof_handler());
      MaterialAveraging::average(cellwise_viscosity_averaging,
                                 current_cell,
                                 fe_values->get_quadrature(),
                                 fe_values->get_mapping(),
                                 MaterialProperties::Property::viscosity,
                                 out);

      // Fill the reaction terms for plastic strain if requested
      if (enable_plastic_strain_weakening && fill_reaction_terms_for_plastic_strain)
        {
          const unsigned int plastic_strain_index = this->introspection().compositional_index_for_name("plastic_strain");
          for (unsigned int i = 0; i < n_evaluation_points; ++i)
            {
              const double plastic_strain_rate = MaterialUtilities::average_value(volume_fractions[i], 
                                                                                  composition_plastic_strain_rates, 
                                                                                  MaterialUtilities::arithmetic);
              out.reaction_terms[i][plastic_strain_index] = plastic_strain_rate * this->get_timestep();
            }
        }
    }



    template <int dim>
    void
    ViscoPlasticSimple<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // Compute the volume fractions, the EoS variables and the average elastic shear modulus
      // (if requested) at each evaluation point.
      std::vector<std::vector<double>> volume_fractions(in.n_evaluation_points());
      EquationOfStateOutputs<dim> eos_outputs(this->introspection().get_number_of_fields_of_type(
        CompositionalFieldDescription::chemical_composition) + 1);
      std::vector<double> average_elastic_shear_moduli(in.n_evaluation_points());

      for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
        {
          volume_fractions[i] = MaterialUtilities::compute_only_composition_fractions(
                                  in.composition[i],
                                  this->introspection().chemical_composition_field_indices());

          equation_of_state.evaluate(in, i, eos_outputs);
          
          out.densities[i] = MaterialUtilities::average_value(volume_fractions[i], 
                                                              eos_outputs.densities, 
                                                              MaterialUtilities::arithmetic);
          out.thermal_expansion_coefficients[i] = MaterialUtilities::average_value(volume_fractions[i], 
                                                                                   eos_outputs.thermal_expansion_coefficients, 
                                                                                   MaterialUtilities::arithmetic);
          out.specific_heat[i] = MaterialUtilities::average_value(volume_fractions[i], 
                                                                  eos_outputs.specific_heat_capacities, 
                                                                  MaterialUtilities::arithmetic);

          if (define_conductivities == false)
            {
              double thermal_diffusivity = 0.0;
              for (unsigned int j = 0; j < volume_fractions[i].size(); ++j)
                thermal_diffusivity += volume_fractions[i][j] * thermal_diffusivities[j];

              // Thermal conductivity at the given positions. If the temperature equation uses
              // the reference density profile formulation, use the reference density to
              // calculate thermal conductivity. Otherwise, use the real density. If the adiabatic
              // conditions are not yet initialized, the real density will still be used.
              if (this->get_parameters().formulation_temperature_equation ==
                  Parameters<dim>::Formulation::TemperatureEquation::reference_density_profile &&
                  this->get_adiabatic_conditions().is_initialized())
                out.thermal_conductivities[i] = thermal_diffusivity * out.specific_heat[i] *
                                                this->get_adiabatic_conditions().density(in.position[i]);
              else
                out.thermal_conductivities[i] = thermal_diffusivity * out.specific_heat[i] * out.densities[i];
            }
          else
            {
              // Use thermal conductivity values specified in the parameter file, if this
              // option was selected.
              out.thermal_conductivities[i] = MaterialUtilities::average_value(volume_fractions[i], 
                                                                               thermal_conductivities, 
                                                                               MaterialUtilities::arithmetic);
            }

          out.compressibilities[i] = MaterialUtilities::average_value(volume_fractions[i], 
                                                                      eos_outputs.compressibilities, 
                                                                      MaterialUtilities::arithmetic);
          out.entropy_derivative_pressure[i] = MaterialUtilities::average_value(volume_fractions[i], 
                                                                                eos_outputs.entropy_derivative_pressure, 
                                                                                MaterialUtilities::arithmetic);
          out.entropy_derivative_temperature[i] = MaterialUtilities::average_value(volume_fractions[i], 
                                                                                   eos_outputs.entropy_derivative_temperature, 
                                                                                   MaterialUtilities::arithmetic);

          if (this->get_parameters().enable_elasticity)
            average_elastic_shear_moduli[i] = MaterialUtilities::average_value(volume_fractions[i],
                                                                               elastic_rheology.get_elastic_shear_moduli(),
                                                                               pointwise_viscosity_averaging);
          // Initialize the reaction terms
          for (unsigned int j = 0; j < in.composition[i].size(); ++j)
            out.reaction_terms[i][j] = 0.;
        }

      // Compute the effective viscosity if requested.
      // Always compute the viscosity if additional outputs are requested, because the viscosity is needed
      // to compute the elastic force term.
      if (in.requests_property(MaterialProperties::viscosity) || in.requests_property(MaterialProperties::additional_outputs))
        {
          double average_strain_rate    = 0.;
          double average_pressure       = 0.;
          double average_temperature    = 0.;
          double average_plastic_strain = 0.;

          double cell_measure = 0.;

          const bool use_reference_strain_rate = this->simulator_is_past_initialization() == false
                                                 ||
                                                 (this->get_timestep_number() == 0 &&
                                                  this->get_nonlinear_iteration() == 0);

          const unsigned int plastic_strain_index = (enable_plastic_strain_weakening ?
                                                     this->introspection().compositional_index_for_name("plastic_strain") :
                                                     numbers::invalid_unsigned_int);
          
          // Only create the FEValues object the first time we get here
          if (fe_values == nullptr)
            {
              // We do not know whether this function is called by Stokes assembler or
              // advection assembler, so the quadratures stored in Introspection are
              // not applicable here. We need to guess the number of quadrature points
              // in one direction.
              unsigned int n_q_points_1d = 1;
              bool bingo = false;
              while (n_q_points_1d < 10)
                {
                  if (in.n_evaluation_points() == Utilities::fixed_power<dim, unsigned int>(n_q_points_1d))
                    {
                      bingo = true;
                      break;
                    }
                  ++n_q_points_1d;
                }
              AssertThrow(bingo, ExcNotImplemented());

              fe_values = std::make_unique<FEValues<dim>>(this->get_mapping(),
                                                          this->get_fe(),
                                                          QGauss<dim>(n_q_points_1d),
                                                          update_JxW_values);
            }

          fe_values->reinit(in.current_cell);

          for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
            {
              double edot_ii_sqr;
              if (use_reference_strain_rate)
                edot_ii_sqr = ref_strain_rate * ref_strain_rate;
              else if (this->get_parameters().enable_elasticity == false)
                {
                  // Calculate the second moment invariant of the deviagtoric strain rate tensor.
                  edot_ii_sqr = std::max(-second_invariant(deviator(in.strain_rate[i])), 0.);
                }
              else
                {
                  // Replace edot_ii the viscoelastic strain rate, which includes a term that accounts for
                  // the elastic stress arising from the previous time step.
                  SymmetricTensor<2,dim> stress_old;
                  for (unsigned int j = 0; j < SymmetricTensor<2,dim>::n_independent_components; ++j)
                    stress_old.access_raw_entry(j) = in.composition[i][j];

                  const SymmetricTensor<2,dim> viscoelastic_strain_rate =
                    elastic_rheology.calculate_viscoelastic_strain_rate(in.strain_rate[i],
                                                                        stress_old,
                                                                        average_elastic_shear_moduli[i]);
                  
                  // Rheology::Elasticity::calculate_viscoelastic_strain_rate() returns the deviatoric
                  // part of the viscoelastic strain rate, so we do not have to call deviator() here
                  edot_ii_sqr = std::max(-second_invariant(viscoelastic_strain_rate), 0.);
                }

              // Choice of activation volume depends on whether there is an adiabatic temperature
              // gradient used when calculating the viscosity. This allows the same activation volume
              // to be used in incompressible and compressible models.
              const double T = this->simulator_is_past_initialization()
                               ?
                               in.temperature[i] + adiabatic_temperature_gradient_for_viscosity * in.pressure[i]
                               :
                               this->get_adiabatic_conditions().temperature(in.position[i]);

              AssertThrow(T != 0, ExcMessage(
                            "The temperature used in the calculation of the visco-plastic rheology is zero. "
                            "This is not allowed, because this value is used to divide through. It is probably "
                            "being caused by the temperature being zero somewhere in the model. The relevant "
                            "values for debugging are: temperature (" + Utilities::to_string(in.temperature[i]) +
                            "), adiabatic_temperature_gradient_for_viscosity ("
                            + Utilities::to_string(adiabatic_temperature_gradient_for_viscosity) + ") and pressure ("
                            + Utilities::to_string(in.pressure[i]) + ")."));

              // Determine whether to use the adiabatic pressure instead of the full pressure (default)
              // when calculating the viscosity.
              const double p = use_adiabatic_pressure_in_viscosity
                               ?
                               this->get_adiabatic_conditions().pressure(in.position[i]) 
                               :
                               in.pressure[i];

              const double e_p = enable_plastic_strain_weakening ?
                                 in.composition[i][plastic_strain_index] :
                                 0.;

              const double JxW = fe_values->JxW(i);

              average_strain_rate    += JxW * edot_ii_sqr;
              average_pressure       += JxW * p;
              average_temperature    += JxW * T;
              average_plastic_strain += JxW * e_p;

              cell_measure += JxW;
            }

          average_strain_rate = std::max(std::sqrt(average_strain_rate / cell_measure), min_strain_rate);
          average_pressure /= cell_measure;
          average_temperature /= cell_measure;
          average_plastic_strain /= cell_measure;

          // Calculate the average viscosity.
          calculate_average_viscosity(volume_fractions,
                                      average_strain_rate, 
                                      average_pressure, 
                                      average_temperature,
                                      average_plastic_strain,
                                      true, out);

          // Compute viscosity derivatives if they are requested
          if (MaterialModel::MaterialModelDerivatives<dim> *derivatives =
              out.template get_additional_output<MaterialModel::MaterialModelDerivatives<dim>>())
            {
              const double finite_difference_accuracy = 1.e-7;

              MaterialModel::MaterialModelOutputs<dim> out_forward(in.n_evaluation_points(),
                                                                   in.composition[0].size());

              // Compute viscosity derivative w.r.t. strain rate
              const double strain_rate_forward = average_strain_rate * (1. + finite_difference_accuracy);
              calculate_average_viscosity(volume_fractions,
                                          strain_rate_forward,
                                          average_pressure,
                                          average_temperature,
                                          average_plastic_strain,
                                          false, out_forward);

              const double deta_deps = (out_forward.viscosities[0] - out.viscosities[0]) /
                                       (strain_rate_forward * finite_difference_accuracy);

              // Compute viscosity derivative w.r.t. pressure
              const double pressure_forward = average_pressure * (1. + finite_difference_accuracy);
              calculate_average_viscosity(volume_fractions,
                                          average_strain_rate,
                                          pressure_forward,
                                          average_temperature,
                                          average_plastic_strain,
                                          false, out_forward);

              const double deta_dp = (std::abs(average_pressure) > std::numeric_limits<double>::min() 
                                      ?
                                      (out_forward.viscosities[0] - out.viscosities[0]) /
                                      (pressure_forward * finite_difference_accuracy) 
                                      :
                                      0.0);
              for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
                {
                  derivatives->viscosity_derivative_wrt_strain_rate[i] = deta_deps;
                  derivatives->viscosity_derivative_wrt_pressure[i] = deta_dp;
                }
            }
        }

      // Finally, fill the reaction terms and additional outputs of the elastic rheology
      if (this->get_parameters().enable_elasticity)
        {
          elastic_rheology.fill_elastic_outputs(in, average_elastic_shear_moduli, out);
          elastic_rheology.fill_reaction_outputs(in, average_elastic_shear_moduli, out);
        }
    }



    template <int dim>
    bool
    ViscoPlasticSimple<dim>::is_compressible() const
    {
      return true;
    }



    template <int dim>
    void
    ViscoPlasticSimple<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Visco plastic simple");
        {
          // Reference and minimum/maximum values
          prm.declare_entry("Reference strain rate", "1.e-15", Patterns::Double(0.),
                            "Reference strain rate for the first time step. Units: \\si{\\per\\second}.");
          prm.declare_entry("Minimum strain rate", "1.e-20", Patterns::Double(0.),
                            "Stabilizes strain dependent viscosity. Units: \\si{\\per\\second}.");
          prm.declare_entry("Minimum viscosity", "1.e17", Patterns::Double(0.),
                            "Lower cutoff for effective viscosity. Units: \\si{\\pascal\\second}.");
          prm.declare_entry("Maximum viscosity", "1.e28", Patterns::Double(0.),
                            "Upper cutoff for effective viscosity. Units: \\si{\\pascal\\second}.");

          // Thermal conductivity
          prm.declare_entry("Define thermal conductivities", "false", Patterns::Bool(),
                            "Whether to directly define thermal conductivities for each compositional field "
                            "instead of calculating the values through the specified thermal diffusivities, "
                            "densities, and heat capacities.");
          prm.declare_entry("Thermal diffusivities", "0.8e-6", Patterns::List(Patterns::Double(0.)),
                            "List of thermal diffusivities, for background material and compositional fields, "
                            "for a total of N+1 values, where N is the number of all compositional fields or only "
                            "those corresponding to chemical compositions. "
                            "If only one value is given, then all use the same value. "
                            "Units: \\si{\\meter\\squared\\per\\second}.");
          prm.declare_entry("Thermal conductivities", "3.0",
                            Patterns::List(Patterns::Double(0)),
                            "List of thermal conductivities, for background material and compositional fields, "
                            "for a total of N+1 values, where N is the number of all compositional fields or only "
                            "those corresponding to chemical compositions. "
                            "If only one value is given, then all use the same value. "
                            "Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");

          // rheological parameters
          prm.declare_entry("Cellwise viscosity averaging scheme", "harmonic average",
                            Patterns::Selection(MaterialAveraging::get_averaging_operation_names()),
                            "The 'visco plastic simple' model requires the viscosity to be "
                            "averaged in each cell."
                            "\n\n"
                            " Possible choices: " 
                            + MaterialAveraging::get_averaging_operation_names()
                            +
                            "\n\n"
                            "Note that the option 'none' is forbidden in this material model.");
          prm.declare_entry("Pointwise viscosity averaging scheme", "harmonic",
                            Patterns::Selection("arithmetic|harmonic|geometric|maximum composition"),
                            "When more than one compositional field is present at a point "
                            "with different viscosities, we need to come up with an average "
                            "viscosity at that point. Select a weighted harmonic, arithmetic, "
                            "geometric, or maximum composition.");
          prm.declare_entry("Adiabat temperature gradient for viscosity", "0.0", Patterns::Double(0.),
                            "Add an adiabatic temperature gradient to the temperature used in the flow law "
                            "so that the activation volume is consistent with what one would use in a "
                            "earth-like (compressible) model. Default is set so this is off. "
                            "Note that this is a linear approximation of the real adiabatic gradient, which "
                            "is okay for the upper mantle, but is not really accurate for the lower mantle. "
                            "Using a pressure gradient of 32436 Pa/m, then a value of "
                            "0.3 K/km = 0.0003 K/m = 9.24e-09 K/Pa gives an earth-like adiabat."
                            "Units: \\si{\\kelvin\\per\\pascal}.");
          prm.declare_entry("Use adiabatic pressure for viscosity", "false", Patterns::Bool(),
                            "Whether to use the adiabatic pressure instead of the full "
                            "pressure (default) when calculating the viscosity. This may be "
                            "helpful in models where the full pressure has unusually large "
                            "variations, resulting in solver convergence issues. Be aware that "
                            "this setting will change the plastic shear band angle.");
          prm.declare_entry("Allow negative pressure in plasticity", "false", Patterns::Bool(),
                            "Whether to allow negative pressures to be used in the computation "
                            "of plastic yield stresses and viscosities. Setting this parameter "
                            "to true may be advantageous in models without gravity where the "
                            "dynamic stresses are much higher than the lithostatic pressure. "
                            "If false, the minimum pressure in the plasticity formulation will "
                            "be set to zero.");
          prm.declare_entry("Viscous flow law", "composite",
                            Patterns::Selection("diffusion|dislocation|composite"),
                            "Select what type of viscosity law to use between diffusion, "
                            "dislocation and composite options.");

          // Plasticity strain weakening
          prm.declare_entry("Enable plasticity strain weakening", "false", Patterns::Bool(),
                            "Wether to apply strain weakening to cohesion and internal angle of friction "
                            "based on accumulated plastic strain.");
          prm.declare_entry("Start plasticity strain weakening intervals", "0.",
                            Patterns::List(Patterns::Double(0.)),
                            "List of strain weakening interval initial strains "
                            "for the cohesion and friction angle parameters of the "
                            "background material and compositional fields, "
                            "for a total of N+1 values, where N is the number of all compositional fields "
                            "or only those corresponding to chemical compositions. "
                            "If only one value is given, then all use the same value. Units: None.");

          prm.declare_entry("End plasticity strain weakening intervals", "1.",
                            Patterns::List(Patterns::Double (0.)),
                            "List of strain weakening interval final strains "
                            "for the cohesion and friction angle parameters of the "
                            "background material and compositional fields, "
                            "for a total of N+1 values, where N is the number of all compositional fields "
                            "or only those corresponding to chemical compositions. "
                            "If only one value is given, then all use the same value.  Units: None.");

          prm.declare_entry("Cohesion strain weakening factors", "1.",
                            Patterns::List(Patterns::Double (0.)),
                            "List of cohesion strain weakening factors "
                            "for background material and compositional fields, "
                            "for a total of N+1 values, where N is the number of all compositional fields "
                            "or only those corresponding to chemical compositions. "
                            "If only one value is given, then all use the same value.  Units: None.");

          prm.declare_entry("Friction strain weakening factors", "1.",
                            Patterns::List(Patterns::Double (0.)),
                            "List of friction strain weakening factors "
                            "for background material and compositional fields, "
                            "for a total of N+1 values, where N is the number of all compositional fields "
                            "or only those corresponding to chemical compositions. "
                            "If only one value is given, then all use the same value.  Units: None.");

          EquationOfState::MulticomponentIncompressible<dim>::declare_parameters(prm);

          Rheology::DiffusionCreep<dim>::declare_parameters(prm);

          Rheology::DislocationCreep<dim>::declare_parameters(prm);

          Rheology::Elasticity<dim>::declare_parameters(prm);

          Rheology::DruckerPrager<dim>::declare_parameters(prm);

          Rheology::FrictionModels<dim>::declare_parameters(prm);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    ViscoPlasticSimple<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Visco plastic simple");
        {
          // Make options file for parsing maps to double arrays
          std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();
          compositional_field_names.insert(compositional_field_names.begin(), "background");

          std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();
          chemical_field_names.insert(chemical_field_names.begin(), "background");

          Utilities::MapParsing::Options options(chemical_field_names, "");
          options.list_of_allowed_keys = compositional_field_names;

          // Reference and minimum/maximum values
          ref_strain_rate = prm.get_double("Reference strain rate");
          min_strain_rate = prm.get_double("Minimum strain rate");
          AssertThrow(ref_strain_rate >= min_strain_rate,
                      ExcMessage("The reference strain rate should be larger than or equal to the minimum strain rate."));

          min_viscosity = prm.get_double("Minimum viscosity");
          max_viscosity = prm.get_double("Maximum viscosity");
          AssertThrow(max_viscosity >= min_viscosity,
                      ExcMessage("The maximum viscosity should be larger than or equal to the minimum viscosity."));
          
          // Thermal conductivity
          define_conductivities = prm.get_bool("Define thermal conductivities");

          options.property_name = "Thermal diffusivities";
          thermal_diffusivities = Utilities::MapParsing::parse_map_to_double_array(prm.get("Thermal diffusivities"), options);

          options.property_name = "Thermal conductivities";
          thermal_conductivities = Utilities::MapParsing::parse_map_to_double_array(prm.get("Thermal conductivities"), options);

          // Rheological parameters
          pointwise_viscosity_averaging = MaterialUtilities::parse_compositional_averaging_operation("Pointwise viscosity averaging scheme", prm);
          cellwise_viscosity_averaging = MaterialAveraging::parse_averaging_operation_name(prm.get("Cellwise viscosity averaging scheme"));
          AssertThrow(cellwise_viscosity_averaging != MaterialAveraging::none,
                      ExcMessage("The material model 'visco plastic simple' requires the viscosity to be averaged in each cell. "
                                 "Please select an averaging scheme for parameter entry 'Cellwise viscosity averaging scheme' "
                                 "other than 'none'."));
          AssertThrow(this->get_parameters().material_averaging == MaterialAveraging::none,
                      ExcMessage("Material model 'visco plastic simple' only works when the parameter entry 'Material averaging' "
                                 "is set to 'none'."));

          use_adiabatic_pressure_in_viscosity = prm.get_bool("Use adiabatic pressure in viscosity");
          allow_negative_pressure_in_plasticity = prm.get_bool("Allow negative pressure in plasticity");
          adiabatic_temperature_gradient_for_viscosity = prm.get_double("Adiabat temperature gradient for viscosity");
          if (this->get_heating_model_manager().adiabatic_heating_enabled())
            AssertThrow (adiabatic_temperature_gradient_for_viscosity == 0.0,
                         ExcMessage("If adiabatic heating is enabled you should not add another adiabatic gradient"
                                    "to the temperature for computing the viscosity, because the ambient"
                                    "temperature profile already includes the adiabatic gradient."));

          if (prm.get("Viscous flow law") == "composite")
            viscous_flow_law = composite;
          else if (prm.get("Viscous flow law") == "diffusion")
            viscous_flow_law = diffusion;
          else if (prm.get("Viscous flow law") == "dislocation")
            viscous_flow_law = dislocation;
          else
            AssertThrow(false, ExcMessage("Not a valid viscous flow law"));

          // Plasticity strain weakening
          enable_plastic_strain_weakening = prm.get_bool("Enable plasticity strain weakening");
          if (enable_plastic_strain_weakening)
            AssertThrow(this->introspection().compositional_name_exists("plastic_strain"),
                        ExcMessage("Plastic strain weakening mechanism only works if there is a "
                                   "compositional field called plastic_strain."));

          options.property_name = "Start plasticity strain weakening intervals";
          start_plastic_strain_weakening_intervals = 
            Utilities::MapParsing::parse_map_to_double_array(prm.get("Start plasticity strain weakening intervals"), options);

          options.property_name = "End plasticity strain weakening intervals";
          end_plastic_strain_weakening_intervals =
            Utilities::MapParsing::parse_map_to_double_array(prm.get("End plasticity strain weakening intervals"), options);

          options.property_name = "Cohesion strain weakening factors";
          cohesion_strain_weakening_factors =
            Utilities::MapParsing::parse_map_to_double_array(prm.get("Cohesion strain weakening factors"), options);

          options.property_name = "Friction strain weakening factors";
          friction_strain_weakening_factors = 
            Utilities::MapParsing::parse_map_to_double_array(prm.get("Friction strain weakening factors"), options);

          equation_of_state.initialize_simulator(this->get_simulator());
          equation_of_state.parse_parameters(prm);

          diffusion_creep.initialize_simulator(this->get_simulator());
          diffusion_creep.parse_parameters(prm);

          dislocation_creep.initialize_simulator(this->get_simulator());
          dislocation_creep.parse_parameters(prm);

          if (this->get_parameters().enable_elasticity)
            {
              elastic_rheology.initialize_simulator(this->get_simulator());
              elastic_rheology.parse_parameters(prm);
            }

          drucker_prager_plasticity.initialize_simulator(this->get_simulator());
          drucker_prager_plasticity.parse_parameters(prm);

          friction_models.initialize_simulator(this->get_simulator());
          friction_models.parse_parameters(prm);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity            = NonlinearDependence::temperature | 
                                                    NonlinearDependence::pressure | 
                                                    NonlinearDependence::strain_rate | 
                                                    NonlinearDependence::compositional_fields;
      this->model_dependence.density              = NonlinearDependence::temperature |
                                                    NonlinearDependence::compositional_fields;
      this->model_dependence.compressibility      = NonlinearDependence::none;
      this->model_dependence.specific_heat        = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::temperature | 
                                                    NonlinearDependence::pressure |
                                                    NonlinearDependence::compositional_fields;
    }



    template <int dim>
    void
    ViscoPlasticSimple<dim>::
    create_additional_named_outputs(MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (this->get_parameters().enable_elasticity)
        elastic_rheology.create_elastic_outputs(out);
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(ViscoPlasticSimple,
                                   "visco plastic simple",
                                   "")
  }
}
