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

#include <aspect/material_model/visco_plastic.h>
#include <aspect/utilities.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/signaling_nan.h>
#include <aspect/newton.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    ViscoPlastic<dim>::initialize()
    {
      if (use_dominant_phase_for_viscosity)
        {
          phase_function_discrete->initialize();
        }
    }


    template <int dim>
    bool
    ViscoPlastic<dim>::
    is_yielding(const MaterialModelInputs<dim> &in) const
    {
      Assert(in.n_evaluation_points() == 1, ExcInternalError());

      const std::vector<double> volume_fractions = MaterialUtilities::compute_only_composition_fractions(in.composition[0],
                                                   this->introspection().chemical_composition_field_indices());

      /* The following handles phases in a similar way as in the 'evaluate' function.
       * Results then enter the calculation of plastic yielding.
       */
      std::vector<double> phase_function_values(phase_function.n_phase_transitions(), 0.0);

      if (phase_function.n_phase_transitions() > 0)
        {
          const double gravity_norm = this->get_gravity_model().gravity_vector(in.position[0]).norm();

          double reference_density;
          if (this->get_adiabatic_conditions().is_initialized())
            {
              reference_density = this->get_adiabatic_conditions().density(in.position[0]);
            }
          else
            {
              EquationOfStateOutputs<dim> eos_outputs_all_phases (n_phases);
              equation_of_state.evaluate(in, 0, eos_outputs_all_phases);
              reference_density = eos_outputs_all_phases.densities[0];
            }

          MaterialUtilities::PhaseFunctionInputs<dim> phase_inputs(in.temperature[0],
                                                                   in.pressure[0],
                                                                   this->get_geometry_model().depth(in.position[0]),
                                                                   gravity_norm*reference_density,
                                                                   numbers::invalid_unsigned_int);

          for (unsigned int j=0; j < phase_function.n_phase_transitions(); ++j)
            {
              phase_inputs.phase_transition_index = j;
              phase_function_values[j] = phase_function.compute_value(phase_inputs);
            }
        }

      /* The following returns whether or not the material is plastically yielding
       * as documented in evaluate.
       */
      const IsostrainViscosities isostrain_viscosities = rheology->calculate_isostrain_viscosities(in, 0, volume_fractions, phase_function_values, phase_function.n_phase_transitions_for_each_composition());

      std::vector<double>::const_iterator max_composition = std::max_element(volume_fractions.begin(), volume_fractions.end());
      const bool plastic_yielding = isostrain_viscosities.composition_yielding[std::distance(volume_fractions.begin(), max_composition)];

      return plastic_yielding;
    }



    template <int dim>
    void
    ViscoPlastic<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      EquationOfStateOutputs<dim> eos_outputs (this->introspection().get_number_of_fields_of_type(CompositionalFieldDescription::chemical_composition)+1);
      EquationOfStateOutputs<dim> eos_outputs_all_phases (n_phases);

      std::vector<double> average_elastic_shear_moduli (in.n_evaluation_points());

      // Store value of phase function for each phase and composition
      // While the number of phases is fixed, the value of the phase function is updated for every point
      std::vector<double> phase_function_values(phase_function.n_phase_transitions(), 0.0);

      std::vector<double> phase_function_discrete_values = (use_dominant_phase_for_viscosity?
                                                            std::vector<double>(phase_function_discrete->n_phase_transitions(), 0.0): std::vector<double>());


      // Loop through all requested points
      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          // First compute the equation of state variables and thermodynamic properties
          equation_of_state.evaluate(in, i, eos_outputs_all_phases);

          const double gravity_norm = this->get_gravity_model().gravity_vector(in.position[i]).norm();
          const double reference_density = (this->get_adiabatic_conditions().is_initialized())
                                           ?
                                           this->get_adiabatic_conditions().density(in.position[i])
                                           :
                                           eos_outputs_all_phases.densities[0];

          // The phase index is set to invalid_unsigned_int, because it is only used internally
          // in phase_average_equation_of_state_outputs to loop over all existing phases
          MaterialUtilities::PhaseFunctionInputs<dim> phase_inputs(in.temperature[i],
                                                                   in.pressure[i],
                                                                   this->get_geometry_model().depth(in.position[i]),
                                                                   gravity_norm*reference_density,
                                                                   numbers::invalid_unsigned_int);

          // Compute value of phase functions
          for (unsigned int j=0; j < phase_function.n_phase_transitions(); ++j)
            {
              phase_inputs.phase_transition_index = j;
              phase_function_values[j] = phase_function.compute_value(phase_inputs);
            }

          // Average by value of gamma function to get value of compositions
          phase_average_equation_of_state_outputs(eos_outputs_all_phases,
                                                  phase_function_values,
                                                  n_phase_transitions_for_each_chemical_composition,
                                                  eos_outputs);

          const std::vector<double> volume_fractions = MaterialUtilities::compute_only_composition_fractions(in.composition[i],
                                                       this->introspection().chemical_composition_field_indices());

          // not strictly correct if thermal expansivities are different, since we are interpreting
          // these compositions as volume fractions, but the error introduced should not be too bad.
          out.densities[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.densities, MaterialUtilities::arithmetic);
          out.thermal_expansion_coefficients[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.thermal_expansion_coefficients, MaterialUtilities::arithmetic);
          out.specific_heat[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.specific_heat_capacities, MaterialUtilities::arithmetic);

          if (define_conductivities == false)
            {
              double thermal_diffusivity = 0.0;

              for (unsigned int j=0; j < volume_fractions.size(); ++j)
                thermal_diffusivity += volume_fractions[j] * thermal_diffusivities[j];

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
              out.thermal_conductivities[i] = MaterialUtilities::average_value (volume_fractions, thermal_conductivities, MaterialUtilities::arithmetic);
            }

          out.compressibilities[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.compressibilities, MaterialUtilities::arithmetic);
          out.entropy_derivative_pressure[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.entropy_derivative_pressure, MaterialUtilities::arithmetic);
          out.entropy_derivative_temperature[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.entropy_derivative_temperature, MaterialUtilities::arithmetic);

          // Compute the effective viscosity if requested and retrieve whether the material is plastically yielding.
          // Also always compute the viscosity if additional outputs are requested, because the viscosity is needed
          // to compute the elastic force term and the elastic reaction rates.
          bool plastic_yielding = false;
          IsostrainViscosities isostrain_viscosities;

          if (in.requests_property(MaterialProperties::viscosity) || in.requests_property(MaterialProperties::additional_outputs) ||
              (this->get_parameters().enable_elasticity && in.requests_property(MaterialProperties::reaction_rates) ))
            {
              // Currently, the viscosities for each of the compositional fields are calculated assuming
              // isostrain amongst all compositions, allowing calculation of the viscosity ratio.
              // TODO: This is only consistent with viscosity averaging if the arithmetic averaging
              // scheme is chosen. It would be useful to have a function to calculate isostress viscosities.
              // Methods for averaging among phases depend on whether the most dominant phases are looked up
              // for each composition in its respective lookup table. A separate phase function class
              // manages the averaging if the user chooses this option.
              if (use_dominant_phase_for_viscosity)
                {
                  for (unsigned int j=0; j < phase_function_discrete->n_phase_transitions(); ++j)
                    {
                      phase_inputs.phase_transition_index = j;
                      phase_function_discrete_values[j] = phase_function_discrete->compute_value(phase_inputs);
                    }
                  isostrain_viscosities =
                    rheology->calculate_isostrain_viscosities(in, i, volume_fractions, phase_function_discrete_values,  phase_function_discrete->n_phase_transitions_for_each_chemical_composition());
                }
              else
                {
                  isostrain_viscosities =
                    rheology->calculate_isostrain_viscosities(in, i, volume_fractions, phase_function_values, n_phase_transitions_for_each_chemical_composition);
                }

              // The isostrain condition implies that the viscosity averaging should be arithmetic (see above).
              // We have given the user freedom to apply alternative bounds, because in diffusion-dominated
              // creep (where n_diff=1) viscosities are stress and strain-rate independent, so the calculation
              // of compositional field viscosities is consistent with any averaging scheme.
              out.viscosities[i] = MaterialUtilities::average_value(volume_fractions, isostrain_viscosities.composition_viscosities, rheology->viscosity_averaging);

              // Decide based on the maximum composition if material is yielding.
              // This avoids for example division by zero for harmonic averaging (as plastic_yielding
              // holds values that are either 0 or 1), but might not be consistent with the viscosity
              // averaging chosen.
              std::vector<double>::const_iterator max_composition = std::max_element(volume_fractions.begin(), volume_fractions.end());
              plastic_yielding = isostrain_viscosities.composition_yielding[std::distance(volume_fractions.begin(), max_composition)];

              // Compute viscosity derivatives if they are requested
              if (const std::shared_ptr<MaterialModel::MaterialModelDerivatives<dim>> derivatives =
                    out.template get_additional_output_object<MaterialModel::MaterialModelDerivatives<dim>>())

                rheology->compute_viscosity_derivatives(i, volume_fractions,
                                                        isostrain_viscosities,
                                                        in, out, phase_function_values,
                                                        n_phase_transitions_for_each_chemical_composition);
            }
          else
            {
              // The viscosity was not requested. Poison its value, along with the other
              // quantities we set above and that would otherwise remain uninitialized
              isostrain_viscosities.composition_yielding.clear();
              isostrain_viscosities.composition_viscosities.clear();
              isostrain_viscosities.drucker_prager_parameters.clear();

              out.viscosities[i] = numbers::signaling_nan<double>();

              if (const std::shared_ptr<MaterialModel::MaterialModelDerivatives<dim>> derivatives =
                    out.template get_additional_output_object<MaterialModel::MaterialModelDerivatives<dim>>())
                {
                  derivatives->viscosity_derivative_wrt_strain_rate[i] = numbers::signaling_nan<SymmetricTensor<2,dim>>();
                  derivatives->viscosity_derivative_wrt_pressure[i] = numbers::signaling_nan<double>();
                }
            }

          // Now compute changes in the compositional fields (e.g., the accumulated strain).
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;

          // Calculate changes in strain invariants and update the reaction terms
          // TODO only when requests_property is set to reaction_terms
          rheology->strain_rheology.fill_reaction_outputs(in, i, rheology->min_strain_rate, plastic_yielding, out);

          // Fill plastic outputs if they exist.
          // The values in isostrain_viscosities only make sense when the calculate_isostrain_viscosities function
          // has been called.
          if (in.requests_property(MaterialProperties::additional_outputs))
            rheology->fill_plastic_outputs(i, volume_fractions, plastic_yielding, in, out, isostrain_viscosities);

          if (this->get_parameters().enable_elasticity)
            {
              // Compute average elastic shear modulus
              average_elastic_shear_moduli[i] = MaterialUtilities::average_value(volume_fractions,
                                                                                 rheology->elastic_rheology.get_elastic_shear_moduli(),
                                                                                 rheology->viscosity_averaging);
            }

          if (const std::shared_ptr<PrescribedPlasticDilation<dim>> plastic_dilation =
                out.template get_additional_output_object<PrescribedPlasticDilation<dim>>())
            {
              const double dilation_lhs_term = MaterialUtilities::average_value(volume_fractions,
                                                                                isostrain_viscosities.dilation_lhs_terms,
                                                                                MaterialUtilities::arithmetic);
              const double dilation_rhs_term = MaterialUtilities::average_value(volume_fractions,
                                                                                isostrain_viscosities.dilation_rhs_terms,
                                                                                MaterialUtilities::arithmetic);

              // When plastic yielding occurs (RHS - LHS * p > 0$), the LHS and RHS terms are set to
              // the values calculated by the Drucker Prager model; otherwise, the LHS and RHS terms
              // should cancel out (RHS = LHS * p) so as to satisfy the loading-unloading conditions.
              plastic_dilation->dilation_lhs_term[i] = dilation_lhs_term;
              if (dilation_rhs_term - dilation_lhs_term * in.pressure[i] > 0)
                {
                  for (unsigned int dim_i = 0; dim_i < dim; ++dim_i)
                    plastic_dilation->dilation_rhs_term[dim_i][i] = dilation_rhs_term;
                }
              else
                {
                  for (unsigned int dim_i = 0; dim_i < dim; ++dim_i)
                    plastic_dilation->dilation_rhs_term[dim_i][i] = dilation_lhs_term * in.pressure[i];
                }
            }
        }

      // If we use the full strain tensor, compute the change in the individual tensor components.
      rheology->strain_rheology.compute_finite_strain_reaction_terms(in, out);

      if (this->get_parameters().enable_elasticity)
        {
          // Fill the elastic outputs with the body force term for the RHS, the viscoelastic strain rate
          // and the viscous dissipation.
          rheology->elastic_rheology.fill_elastic_outputs(in, average_elastic_shear_moduli, out);
          // Fill the elastic additional outputs with the shear modulus, elastic viscosity
          // and deviatoric stress of the current timestep.
          // TODO requests_property is already checked in the fill_ function,
          // but we can also do it here
          //if (in.requests_property(MaterialProperties::additional_outputs))
          rheology->elastic_rheology.fill_elastic_additional_outputs(in, average_elastic_shear_moduli, out);
          // Fill the reaction terms that account for the rotation of the stresses.
          rheology->elastic_rheology.fill_reaction_outputs(in, average_elastic_shear_moduli, out);
          // Fill the reaction_rates that apply the stress update of the previous
          // timestep to the advected and rotated stress computed in the previous timestep ($\tau^{0}$)
          // to obtain $\tau^{t}$.
          rheology->elastic_rheology.fill_reaction_rates(in, average_elastic_shear_moduli, out);
        }
    }



    template <int dim>
    bool
    ViscoPlastic<dim>::
    is_compressible () const
    {
      return equation_of_state.is_compressible();
    }



    template <int dim>
    double ViscoPlastic<dim>::
    get_min_strain_rate () const
    {
      return rheology->min_strain_rate;
    }



    template <int dim>
    void
    ViscoPlastic<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Visco Plastic");
        {
          prm.declare_entry ("Use dominant phase for viscosity","false",
                             Patterns::Bool (),
                             "Whether to look up the dominant phase for each composition in its respective "
                             "material data file to calculate viscosity. This allows each phase to have distinct "
                             "rheological parameterizations.");

          MaterialUtilities::PhaseFunction<dim>::declare_parameters(prm);

          MaterialUtilities::PhaseFunctionDiscrete<dim>::declare_parameters(prm);

          EquationOfState::MulticomponentIncompressible<dim>::declare_parameters (prm);

          Rheology::ViscoPlastic<dim>::declare_parameters(prm);

          // Equation of state parameters
          prm.declare_entry ("Thermal diffusivities", "0.8e-6",
                             Patterns::List(Patterns::Double (0.)),
                             "List of thermal diffusivities, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of all compositional fields or only "
                             "those corresponding to chemical compositions. "
                             "If only one value is given, then all use the same value.  "
                             "Units: \\si{\\meter\\squared\\per\\second}.");
          prm.declare_entry ("Define thermal conductivities","false",
                             Patterns::Bool (),
                             "Whether to directly define thermal conductivities for each compositional field "
                             "instead of calculating the values through the specified thermal diffusivities, "
                             "densities, and heat capacities. ");
          prm.declare_entry ("Thermal conductivities", "3.0",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal conductivities, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of all compositional fields or only "
                             "those corresponding to chemical compositions. "
                             "If only one value is given, then all use the same value. "
                             "Units: $\\frac{\\text{W}{\\text{m}\\text{K}}$.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    ViscoPlastic<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Visco Plastic");
        {
          // Phase transition parameters
          phase_function.initialize_simulator (this->get_simulator());
          phase_function.parse_parameters (prm);

          const std::vector<unsigned int> n_phases_for_each_chemical_composition = phase_function.n_phases_for_each_chemical_composition();
          n_phase_transitions_for_each_chemical_composition = phase_function.n_phase_transitions_for_each_chemical_composition();
          n_phases = phase_function.n_phases_over_all_chemical_compositions();

          // Equation of state parameters
          equation_of_state.initialize_simulator (this->get_simulator());
          equation_of_state.parse_parameters (prm,
                                              std::make_unique<std::vector<unsigned int>>(n_phases_for_each_chemical_composition));

          // Make options file for parsing maps to double arrays
          std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();
          chemical_field_names.insert(chemical_field_names.begin(),"background");

          std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();
          compositional_field_names.insert(compositional_field_names.begin(),"background");

          Utilities::MapParsing::Options options(chemical_field_names, "Thermal diffusivities");
          options.list_of_allowed_keys = compositional_field_names;
          options.allow_multiple_values_per_key = true;
          options.n_values_per_key = n_phases_for_each_chemical_composition;
          options.check_values_per_key = (options.n_values_per_key.size() != 0);
          options.store_values_per_key = (options.n_values_per_key.size() == 0);

          thermal_diffusivities = Utilities::MapParsing::parse_map_to_double_array(prm.get("Thermal diffusivities"), options);
          define_conductivities = prm.get_bool ("Define thermal conductivities");

          options.property_name = "Thermal conductivities";
          thermal_conductivities = Utilities::MapParsing::parse_map_to_double_array (prm.get("Thermal conductivities"), options);

          rheology = std::make_unique<Rheology::ViscoPlastic<dim>>();
          rheology->initialize_simulator (this->get_simulator());

          use_dominant_phase_for_viscosity = prm.get_bool ("Use dominant phase for viscosity");
          if (use_dominant_phase_for_viscosity)
            {
              phase_function_discrete = std::make_unique<MaterialUtilities::PhaseFunctionDiscrete<dim>>();
              phase_function_discrete->initialize_simulator (this->get_simulator());
              phase_function_discrete->parse_parameters (prm);
              rheology->parse_parameters(prm, std::make_unique<std::vector<unsigned int>>(phase_function_discrete->n_phases_for_each_chemical_composition()));
            }
          else
            {
              rheology->parse_parameters(prm, std::make_unique<std::vector<unsigned int>>(n_phases_for_each_chemical_composition));
            }

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



    template <int dim>
    void
    ViscoPlastic<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      rheology->create_plastic_outputs(out);

      if (this->get_parameters().enable_elasticity)
        rheology->elastic_rheology.create_elastic_additional_outputs(out);
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
                                   "An implementation of an incompressible visco(elastic)-plastic rheology "
                                   "with options for selecting dislocation creep, diffusion creep or "
                                   "composite viscous flow laws. Prior to yielding, one may select to "
                                   "modify the viscosity to account for viscoelastic effects by setting the "
                                   "parameter 'Enable elasticity' in subsection Formulation to true. Plasticity "
                                   "limits viscous stresses through a Drucker Prager yield criterion. "
                                   "The implementation of this material model is based heavily on the "
                                   "`DiffusionDislocation' (Bob Myhill), `DruckerPrager' "
                                   "(Anne Glerum), and `Viscoelastic' (John Naliboff) material models. "
                                   "\n\n "
                                   "The viscosity for dislocation or diffusion creep is defined as "
                                   "$ \\eta = \\frac 12 A^{-\\frac{1}{n}} d^{\\frac{m}{n}} "
                                   "\\dot{\\varepsilon}_{ii}^{\\frac{1-n}{n}} "
                                   "\\exp\\left(\\frac{E + PV}{nRT}\\right)$ "
                                   "where $A$ is the prefactor, $n$ is the stress exponent, "
                                   "$\\dot{\\varepsilon}_{ii}$ is the square root of the deviatoric "
                                   "strain rate tensor second invariant, $d$ is grain size, "
                                   "$m$ is the grain size exponent, $E$ is activation energy, "
                                   "$V$ is activation volume, $P$ is pressure, $R$ is the gas "
                                   "exponent and $T$ is temperature. "
                                   "This form of the viscosity equation is commonly used in "
                                   "geodynamic simulations. See, for example, Billen and Hirth "
                                   "(2007), G3, 8, Q08012. Significantly, other studies may use "
                                   "slightly different forms of the viscosity equation leading to "
                                   "variations in how specific terms are defined or combined. For "
                                   "example, the grain size exponent should always be positive in "
                                   "the diffusion viscosity equation used here, while other studies "
                                   "place the grain size term in the denominator and invert the sign "
                                   "of the grain size exponent. When examining previous work, one "
                                   "should carefully check how the viscous prefactor and grain size "
                                   "terms are defined. "
                                   "\n\n "
                                   "One may select to use the diffusion ($\\eta_{\\text{diff}}$; $n=1$, $m\\neq 0$), "
                                   "dislocation ($\\eta_{\\text{disl}}$, $n>1$, $m=0$) or composite "
                                   "$\\frac{\\eta_{\\text{diff}} \\eta_{\\text{disl}}}{\\eta_{\\text{diff}}+\\eta_{\\text{disl}}}$ equation form. "
                                   "\n\n "
                                   "The diffusion and dislocation prefactors can be weakened with a factor "
                                   "between 0 and 1 according to the total or the viscous strain only. "
                                   "\n\n "
                                   "Viscosity is limited through one of two different `yielding' mechanisms. "
                                   "\n\n"
                                   "The first plasticity mechanism limits viscous stress through a "
                                   "Drucker Prager yield criterion, where the yield stress in 3d is  "
                                   "$\\sigma_y = \\frac{6C\\cos(\\phi) + 6P\\sin(\\phi)} "
                                   "{\\sqrt{3}(3+\\sin(\\phi))}$ "
                                   "and "
                                   "$\\sigma_y = C\\cos(\\phi) + P\\sin(\\phi)$ "
                                   "in 2d. Above, $C$ is cohesion and $\\phi$  is the angle of "
                                   "internal friction.  Note that the 2d form is equivalent to the "
                                   "Mohr Coulomb yield surface.  If $\\phi$ is 0, the yield stress "
                                   "is fixed and equal to the cohesion (Von Mises yield criterion). "
                                   "When the viscous stress ($2\\eta {\\varepsilon}_{ii}$) exceeds "
                                   "the yield stress, the viscosity is rescaled back to the yield "
                                   "surface: $\\eta_{y}=\\sigma_{y}/(2{\\varepsilon}_{ii})$. "
                                   "This form of plasticity is commonly used in geodynamic models. "
                                   "See, for example, Thieulot, C. (2011), PEPI 188, pp. 47-68. "
                                   "\n\n"
                                   "The user has the option to linearly reduce the cohesion and "
                                   "internal friction angle as a function of the finite strain magnitude. "
                                   "The finite strain invariant or full strain tensor is calculated through "
                                   "compositional fields within the material model. This implementation is "
                                   "identical to the compositional field finite strain plugin and cookbook "
                                   "described in the manual (author: Gassmoeller, Dannberg). If the user selects to track "
                                   "the finite strain invariant ($e_{ii}$), a single compositional field tracks "
                                   "the value derived from $e_{ii}^t = (e_{ii})^{(t-1)} + \\dot{e}_{ii}\\; dt$, where $t$ and $t-1$ "
                                   "are the current and prior time steps, $\\dot{e}_{ii}$ is the second invariant of the "
                                   "strain rate tensor and $dt$ is the time step size. In the case of the "
                                   "full strain tensor $F$, the finite strain magnitude is derived from the "
                                   "second invariant of the symmetric stretching tensor $L$, where "
                                   "$L = F [F]^T$. The user must specify a single compositional "
                                   "field for the finite strain invariant or multiple fields (4 in 2d, 9 in 3d) "
                                   "for the finite strain tensor. These field(s) must be the first listed "
                                   "compositional fields in the parameter file. Note that one or more of the finite strain "
                                   "tensor components must be assigned a non-zero value initially. This value can be "
                                   "be quite small (e.g., 1.e-8), but still non-zero. While the option to track and use "
                                   "the full finite strain tensor exists, tracking the associated compositional fields "
                                   "is computationally expensive in 3d. Similarly, the finite strain magnitudes "
                                   "may in fact decrease if the orientation of the deformation field switches "
                                   "through time. Consequently, the ideal solution is track the finite strain "
                                   "invariant (single compositional) field within the material and track "
                                   "the full finite strain tensor through particles."
                                   "When only the second invariant of the strain is tracked, one has the option to "
                                   "track the full strain or only the plastic strain. In the latter case, strain is only tracked "
                                   "in case the material is plastically yielding, i.e. the viscous stress > yield stress. "
                                   "\n\n"
                                   "Viscous stress may also be limited by a non-linear stress limiter "
                                   "that has a form similar to the Peierls creep mechanism. "
                                   "This stress limiter assigns an effective viscosity "
                                   "$\\sigma_{\\text{eff}} = \\frac{\\tau_y}{2\\varepsilon_y} "
                                   "{\\frac{\\varepsilon_{ii}}{\\varepsilon_y}}^{\\frac{1}{n_y}-1}$ "
                                   "Above $\\tau_y$ is a yield stress, $\\varepsilon_y$ is the "
                                   "reference strain rate, $\\varepsilon_{ii}$ is the strain rate "
                                   "and $n_y$ is the stress limiter exponent.  The yield stress, "
                                   "$\\tau_y$, is defined through the Drucker Prager yield criterion "
                                   "formulation. This method of limiting viscous stress has been used "
                                   "in various forms within the geodynamic literature \\cite{chri92,vavv02,cibi13,cibi15}."
                                   "When $n_y$ is 1, it essentially becomes a linear viscosity model, "
                                   "and in the limit $n_y\\rightarrow \\infty$ it converges to the "
                                   "standard viscosity rescaling method (concretely, values $n_y>20$ "
                                   "are large enough)."
                                   "\n\n "
                                   "The visco-plastic rheology described above may also be modified to include "
                                   "viscoelastic deformation, thus producing a viscoelastic plastic constitutive "
                                   "relationship. "
                                   "\n\n "
                                   "The viscoelastic rheology behavior takes into account the elastic shear "
                                   "strength (e.g., shear modulus), while the tensile and volumetric "
                                   "strength (e.g., Young's and bulk modulus) are not considered. The "
                                   "model is incompressible and allows specifying an arbitrary number "
                                   "of compositional fields, where each field represents a different "
                                   "rock type or component of the viscoelastic stress tensor. The stress "
                                   "tensor in 2d and 3d, respectively, contains 3 or 6 components. The "
                                   "compositional fields representing these components must be named "
                                   "and listed in a very specific format, which is designed to minimize "
                                   "mislabeling stress tensor components as distinct 'compositional "
                                   "rock types' (or vice versa). For 2d models, three plus three consecutive "
                                   "compositional fields must be labeled 'stress\\_xx', 'stress\\_yy', 'stress\\_xy', "
                                   "'stress\\_xx\\_old', 'stress\\_yy\\_old', and 'stress\\_xy\\_old'. "
                                   "In 3d, six plus six compositional fields must be labeled 'stress\\_xx', "
                                   "'stress\\_yy', 'stress\\_zz', 'stress\\_xy', 'stress\\_xz', 'stress\\_yz', "
                                   "'stress\\_xx\\_old', 'stress\\_yy\\_old', 'stress\\_zz\\_old', 'stress\\_xy\\_old', "
                                   "'stress\\_xz\\_old', 'stress\\_yz\\_old'. "
                                   "\n\n "
                                   "Combining this viscoelasticity implementation with non-linear viscous flow "
                                   "and plasticity produces a constitutive relationship commonly referred to "
                                   "as partial elastoviscoplastic (e.g., pEVP) in the geodynamics community. "
                                   "While extensively discussed and applied within the geodynamics "
                                   "literature, notable references include: "
                                   "Moresi et al. (2003), J. Comp. Phys., v. 184, p. 476-497. "
                                   "Gerya and Yuen (2007), Phys. Earth. Planet. Inter., v. 163, p. 83-105. "
                                   "Gerya (2010), Introduction to Numerical Geodynamic Modeling. "
                                   "Kaus (2010), Tectonophysics, v. 484, p. 36-47. "
                                   "Choi et al. (2013), J. Geophys. Res., v. 118, p. 2429-2444. "
                                   "Keller et al. (2013), Geophys. J. Int., v. 195, p. 1406-1442. "
                                   "\n\n "
                                   "The overview below directly follows Moresi et al. (2003) eqns. 23-38. "
                                   "However, an important distinction between this material model and "
                                   "the studies above is the option to use compositional fields, rather than "
                                   "particles, to track individual components of the viscoelastic stress "
                                   "tensor. Calculating viscoelastic stresses with particles is also implemented, "
                                   "and can be switched on by using particles with the particle property 'elastic stress'. "
                                   "\n\n "
                                   "Moresi et al. (2003) begins (eqn. 23) by writing the deviatoric "
                                   "rate of deformation ($\\hat{D}$) as the sum of elastic "
                                   "($\\hat{D_{e}}$) and viscous ($\\hat{D_{v}}$) components: "
                                   "$\\hat{D} = \\hat{D_{e}} + \\hat{D_{v}}$.  "
                                   "These terms further decompose into "
                                   "$\\hat{D_{v}} = \\frac{\\tau}{2\\eta}$ and "
                                   "$\\hat{D_{e}} = \\frac{\\overset{\\nabla}{\\tau}}{2\\mu}$, where "
                                   "$\\tau$ is the viscous deviatoric stress, $\\eta$ is the shear viscosity, "
                                   "$\\mu$ is the shear modulus and $\\overset{\\nabla}{\\tau}$ is the "
                                   "Jaumann corotational stress rate. This later term (eqn. 24) contains the "
                                   "time derivative of the deviatoric stress ($\\dot{\\tau}$) and terms that "
                                   "account for material spin (e.g., rotation) due to advection: "
                                   "$\\overset{\\nabla}{\\tau} = \\dot{\\tau} + {\\tau}W -W\\tau$. "
                                   "Above, $W$ is the material spin tensor (eqn. 25): "
                                   "$W_{ij} = \\frac{1}{2} \\left (\\frac{\\partial V_{i}}{\\partial x_{j}} - "
                                   "\\frac{\\partial V_{j}}{\\partial x_{i}} \\right )$. "
                                   "\n\n "
                                   "If plasticity is included, the deviatoric rate of deformation may be written as: "
                                   "$\\hat{D} = \\hat{D_{e}} + \\hat{D_{v}} + \\hat{D_{p}}$, where $\\hat{D_{p}}$ "
                                   "is the plastic component. $\\hat{D_{p}}$ decomposes to $\\frac{\\tau_{y}}{2\\eta_{y}}$, "
                                   "where $\\tau_{y}$ is the yield stress and $\\eta_{y}$ is the viscosity rescaled "
                                   "to the yield surface. "
                                   "The Jaumann stress-rate can also be approximated using terms from the "
                                   "previous time step ($t$) and current time step ($t + \\Delta t^{e}$): "
                                   "$\\smash[t]{\\overset{\\nabla}{\\tau}}^{t + \\Delta t^{e}} \\approx "
                                   "\\frac{\\tau^{t + \\Delta t^{e} - \\tau^{t}}}{\\Delta t^{e}} - "
                                   "W^{t}\\tau^{t} + \\tau^{t}W^{t}$. "
                                   "In this material model, the size of the time step above ($\\Delta t^{e}$) "
                                   "can be specified as the numerical time step size or an independent fixed time "
                                   "step. If the latter case is selected, the user has an option to apply a "
                                   "stress averaging scheme to account for the differences between the numerical "
                                   "and fixed elastic time step (eqn. 32). If one selects to use a fixed elastic time "
                                   "step throughout the model run, this can still be achieved by using CFL and "
                                   "maximum time step values that restrict the numerical time step to a specific time."
                                   "\n\n "
                                   "The formulation above allows rewriting the total rate of deformation (eqn. 29) as\n "
                                   "$\\tau^{t + \\Delta t^{e}} = \\eta_{eff} \\left ( "
                                   "2\\hat{D}^{t + \\triangle t^{e}} + \\frac{\\tau^{t}}{\\mu \\Delta t^{e}} + "
                                   "\\frac{W^{t}\\tau^{t} - \\tau^{t}W^{t}}{\\mu}  \\right )$. "
                                   "\n\n "
                                   "The effective viscosity (eqn. 28) is a function of the viscosity ($\\eta$), "
                                   "elastic time step size ($\\Delta t^{e}$) and shear relaxation time "
                                   "($ \\alpha = \\frac{\\eta}{\\mu} $): "
                                   "$\\eta_{eff} = \\eta \\frac{\\Delta t^{e}}{\\Delta t^{e} + \\alpha}$ "
                                   "The magnitude of the shear modulus thus controls how much the effective "
                                   "viscosity is reduced relative to the initial viscosity. "
                                   "\n\n "
                                   "Elastic effects are introduced into the governing Stokes equations through "
                                   "an elastic force term (eqn. 30 updated to the term in eqn. 5 in Farrington et al. 2014) "
                                   "using stresses from the previous time step rotated and advected into the current time step: "
                                   "$F^{e,t} = -\\frac{\\eta_{eff}}{\\mu \\Delta t^{e}} \\tau^{0adv}$. "
                                   "This force term is added onto the right-hand side force vector in the "
                                   "system of equations. "
                                   "\n\n "
                                   "When plastic yielding occurs, the effective viscosity in equation 29 and 30 is the "
                                   "plastic viscosity (equation 36). If the current stress is below the plastic "
                                   "yield stress, the effective viscosity is still as defined in equation 28. "
                                   "During non-linear iterations, we define the current stress prior to yielding "
                                   "(e.g., value compared to yield stress) as "
                                   "$\\tau^{t + \\Delta t^{e}} = \\eta_{eff} \\left ( 2\\hat{D}^{t + \\triangle t^{e}} + "
                                   "\\frac{\\tau^{t}}{\\mu \\Delta t^{e}} \\right ) $"
                                   "\n\n "
                                   "Compositional fields can each be assigned individual values of "
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
