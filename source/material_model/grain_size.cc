/*
  Copyright (C) 2014 - 2024 by the authors of the ASPECT code.

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


#include <aspect/material_model/grain_size.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/material_model/rheology/visco_plastic.h>
#include <aspect/utilities.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace
    {
      std::vector<std::string> make_dislocation_viscosity_outputs_names()
      {
        std::vector<std::string> names;
        names.emplace_back("dislocation_viscosity");
        names.emplace_back("diffusion_viscosity");
        return names;
      }
    }



    template <int dim>
    DislocationViscosityOutputs<dim>::DislocationViscosityOutputs (const unsigned int n_points)
      :
      NamedAdditionalMaterialOutputs<dim>(make_dislocation_viscosity_outputs_names()),
      dislocation_viscosities(n_points, numbers::signaling_nan<double>()),
      diffusion_viscosities(n_points, numbers::signaling_nan<double>())
    {}



    template <int dim>
    std::vector<double>
    DislocationViscosityOutputs<dim>::get_nth_output(const unsigned int idx) const
    {
      AssertIndexRange (idx, 2);
      switch (idx)
        {
          case 0:
            return dislocation_viscosities;

          case 1:
            return diffusion_viscosities;

          default:
            AssertThrow(false, ExcInternalError());
        }
      // we will never get here, so just return something
      return dislocation_viscosities;
    }



    template <int dim>
    void
    GrainSize<dim>::initialize()
    {
      CitationInfo::add("grainsize");

      n_material_data = material_file_names.size();
      for (unsigned i = 0; i < n_material_data; ++i)
        {
          if (material_file_format == perplex)
            material_lookup
            .push_back(std::make_unique<MaterialModel::MaterialUtilities::Lookup::PerplexReader>(datadirectory+material_file_names[i],
                       use_bilinear_interpolation,
                       this->get_mpi_communicator()));
          else if (material_file_format == hefesto)
            material_lookup
            .push_back(std::make_unique<MaterialModel::MaterialUtilities::Lookup::HeFESToReader>(datadirectory+material_file_names[i],
                       datadirectory+derivatives_file_names[i],
                       use_bilinear_interpolation,
                       this->get_mpi_communicator()));
          else
            AssertThrow (false, ExcNotImplemented());
        }
    }



    template <int dim>
    unsigned int
    GrainSize<dim>::
    get_phase_index (const MaterialUtilities::PhaseFunctionInputs<dim> &in) const
    {
      // Since phase transition depth increases monotonically, we only need
      // to check for the first phase that has not yet undergone the transition
      // (phase function value lower than 0.5).
      for (unsigned int j=0; j<n_phase_transitions; ++j)
        {
          MaterialUtilities::PhaseFunctionInputs<dim> phase_inputs(in.temperature,
                                                                   in.pressure,
                                                                   in.depth,
                                                                   in.pressure_depth_derivative,
                                                                   j);

          if (phase_function->compute_value(phase_inputs) < 0.5)
            return j;
        }

      return n_phase_transitions;
    }



    template <int dim>
    double
    GrainSize<dim>::
    diffusion_viscosity (const double temperature,
                         const double adiabatic_temperature,
                         const double adiabatic_pressure,
                         const double grain_size,
                         const double second_strain_rate_invariant,
                         const unsigned int phase_index) const
    {
      double energy_term = std::exp((diffusion_activation_energy[phase_index] + diffusion_activation_volume[phase_index] * adiabatic_pressure)
                                    / (diffusion_creep_exponent[phase_index] * constants::gas_constant * temperature));

      // If the adiabatic profile is already calculated we can use it to limit
      // variations in viscosity due to temperature.
      if (this->get_adiabatic_conditions().is_initialized())
        {
          const double adiabatic_energy_term
            = std::exp((diffusion_activation_energy[phase_index] + diffusion_activation_volume[phase_index] * adiabatic_pressure)
                       / (diffusion_creep_exponent[phase_index] * constants::gas_constant * adiabatic_temperature));

          const double temperature_dependence = energy_term / adiabatic_energy_term;
          if (temperature_dependence > max_temperature_dependence_of_eta)
            energy_term = adiabatic_energy_term * max_temperature_dependence_of_eta;
          if (temperature_dependence < 1.0 / max_temperature_dependence_of_eta)
            energy_term = adiabatic_energy_term / max_temperature_dependence_of_eta;
        }

      const double strain_rate_dependence = (1.0 - diffusion_creep_exponent[phase_index]) / diffusion_creep_exponent[phase_index];

      return diffusion_creep_prefactor[phase_index]
             * std::pow(second_strain_rate_invariant,strain_rate_dependence)
             * std::pow(grain_size, diffusion_creep_grain_size_exponent[phase_index]/diffusion_creep_exponent[phase_index])
             * energy_term;
    }



    template <int dim>
    double
    GrainSize<dim>::
    dislocation_viscosity (const double temperature,
                           const double adiabatic_temperature,
                           const double adiabatic_pressure,
                           const SymmetricTensor<2,dim> &strain_rate,
                           const unsigned int phase_index,
                           const double diffusion_viscosity,
                           const double viscosity_guess) const
    {
      // find out in which phase we are
      double energy_term = std::exp((dislocation_activation_energy[phase_index] + dislocation_activation_volume[phase_index] * adiabatic_pressure)
                                    / (dislocation_creep_exponent[phase_index] * constants::gas_constant * temperature));

      // If we are past the initialization of the adiabatic profile, use it to
      // limit viscosity variations due to temperature.
      if (this->get_adiabatic_conditions().is_initialized())
        {
          const double adiabatic_energy_term
            = std::exp((dislocation_activation_energy[phase_index] + dislocation_activation_volume[phase_index] * adiabatic_pressure)
                       / (dislocation_creep_exponent[phase_index] * constants::gas_constant * adiabatic_temperature));

          const double temperature_dependence = energy_term / adiabatic_energy_term;
          if (temperature_dependence > max_temperature_dependence_of_eta)
            energy_term = adiabatic_energy_term * max_temperature_dependence_of_eta;
          if (temperature_dependence < 1.0 / max_temperature_dependence_of_eta)
            energy_term = adiabatic_energy_term / max_temperature_dependence_of_eta;
        }

      const double strain_rate_dependence = (1.0 - dislocation_creep_exponent[phase_index]) / dislocation_creep_exponent[phase_index];
      const SymmetricTensor<2,dim> shear_strain_rate = strain_rate - 1./dim * trace(strain_rate) * unit_symmetric_tensor<dim>();
      const double second_strain_rate_invariant = std::sqrt(std::max(-second_invariant(shear_strain_rate), 0.));

      // If the strain rate is zero, the dislocation viscosity is infinity.
      if (second_strain_rate_invariant <= std::numeric_limits<double>::min())
        return std::numeric_limits<double>::min();

      // Start the iteration with the full strain rate
      double dis_viscosity;
      if (viscosity_guess == 0)
        dis_viscosity = dislocation_creep_prefactor[phase_index]
                        * std::pow(second_strain_rate_invariant,strain_rate_dependence)
                        * energy_term;
      else
        dis_viscosity = viscosity_guess;

      double dis_viscosity_old = 0;
      unsigned int i = 0;
      while ((std::abs((dis_viscosity-dis_viscosity_old) / dis_viscosity) > dislocation_viscosity_iteration_threshold)
             && (i < dislocation_viscosity_iteration_number))
        {
          const SymmetricTensor<2,dim> dislocation_strain_rate = diffusion_viscosity
                                                                 / (diffusion_viscosity + dis_viscosity) * shear_strain_rate;
          const double dislocation_strain_rate_invariant = std::sqrt(std::max(-second_invariant(dislocation_strain_rate), 0.));

          dis_viscosity_old = dis_viscosity;
          dis_viscosity = dislocation_creep_prefactor[phase_index]
                          * std::pow(dislocation_strain_rate_invariant,strain_rate_dependence)
                          * energy_term;
          ++i;
        }

      Assert(i<dislocation_viscosity_iteration_number,ExcInternalError());

      return dis_viscosity;
    }



    template <int dim>
    double
    GrainSize<dim>::
    enthalpy (const double      temperature,
              const double      pressure,
              const std::vector<double> &compositional_fields,
              const Point<dim> &/*position*/) const
    {
      AssertThrow ((reference_compressibility != 0.0) || use_table_properties,
                   ExcMessage("Currently only compressible models are supported."));

      double enthalpy = 0.0;
      if (n_material_data == 1)
        enthalpy = material_lookup[0]->enthalpy(temperature,pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; ++i)
            enthalpy += compositional_fields[i] * material_lookup[i]->enthalpy(temperature,pressure);
        }
      return enthalpy;
    }



    template <int dim>
    double
    GrainSize<dim>::
    seismic_Vp (const double      temperature,
                const double      pressure,
                const std::vector<double> &compositional_fields,
                const Point<dim> &/*position*/) const
    {
      AssertThrow ((reference_compressibility != 0.0) || use_table_properties,
                   ExcMessage("Currently only compressible models are supported."));

      double vp = 0.0;
      if (n_material_data == 1)
        vp = material_lookup[0]->seismic_Vp(temperature,pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; ++i)
            vp += compositional_fields[i] * material_lookup[i]->seismic_Vp(temperature,pressure);
        }
      return vp;
    }



    template <int dim>
    double
    GrainSize<dim>::
    seismic_Vs (const double      temperature,
                const double      pressure,
                const std::vector<double> &compositional_fields,
                const Point<dim> &/*position*/) const
    {
      AssertThrow ((reference_compressibility != 0.0) || use_table_properties,
                   ExcMessage("Currently only compressible models are supported."));

      double vs = 0.0;
      if (n_material_data == 1)
        vs = material_lookup[0]->seismic_Vs(temperature,pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; ++i)
            vs += compositional_fields[i] * material_lookup[i]->seismic_Vs(temperature,pressure);
        }
      return vs;
    }



    template <int dim>
    double
    GrainSize<dim>::
    density (const double temperature,
             const double pressure,
             const std::vector<double> &compositional_fields,
             const Point<dim> &/*position*/) const
    {
      if (!use_table_properties)
        {
          return reference_rho * std::exp(reference_compressibility * (pressure - this->get_surface_pressure()))
                 * (1 - thermal_alpha * (temperature - reference_T));
        }
      else
        {
          double rho = 0.0;
          if (n_material_data == 1)
            {
              rho = material_lookup[0]->density(temperature,pressure);
            }
          else
            {
              for (unsigned i = 0; i < n_material_data; ++i)
                rho += compositional_fields[i] * material_lookup[i]->density(temperature,pressure);
            }

          return rho;
        }
    }



    template <int dim>
    bool
    GrainSize<dim>::
    is_compressible () const
    {
      return (reference_compressibility != 0)
             || use_table_properties;
    }



    template <int dim>
    double
    GrainSize<dim>::
    compressibility (const double temperature,
                     const double pressure,
                     const std::vector<double> &compositional_fields,
                     const Point<dim> &position) const
    {
      if (!use_table_properties)
        return reference_compressibility;

      double dRhodp = 0.0;
      if (n_material_data == 1)
        dRhodp = material_lookup[0]->dRhodp(temperature,pressure);
      else
        {
          for (unsigned i = 0; i < n_material_data; ++i)
            dRhodp += compositional_fields[i] * material_lookup[i]->dRhodp(temperature,pressure);
        }
      const double rho = density(temperature,pressure,compositional_fields,position);
      return (1/rho)*dRhodp;
    }



    template <int dim>
    double
    GrainSize<dim>::
    thermal_expansion_coefficient (const double      temperature,
                                   const double      pressure,
                                   const std::vector<double> &compositional_fields,
                                   const Point<dim> &/*position*/) const
    {
      double alpha = 0.0;
      if (!use_table_properties)
        return thermal_alpha;
      else
        {
          if (n_material_data == 1)
            alpha = material_lookup[0]->thermal_expansivity(temperature,pressure);
          else
            {
              for (unsigned i = 0; i < n_material_data; ++i)
                alpha += compositional_fields[i] * material_lookup[i]->thermal_expansivity(temperature,pressure);
            }
        }
      alpha = std::max(std::min(alpha,max_thermal_expansivity),min_thermal_expansivity);
      return alpha;
    }



    template <int dim>
    double
    GrainSize<dim>::
    specific_heat (const double temperature,
                   const double pressure,
                   const std::vector<double> &compositional_fields,
                   const Point<dim> &/*position*/) const
    {
      double cp = 0.0;
      if (!use_table_properties)
        return reference_specific_heat;
      else
        {
          if (n_material_data == 1)
            cp = material_lookup[0]->specific_heat(temperature,pressure);
          else
            {
              for (unsigned i = 0; i < n_material_data; ++i)
                cp += compositional_fields[i] * material_lookup[i]->specific_heat(temperature,pressure);
            }
        }
      cp = std::max(std::min(cp,max_specific_heat),min_specific_heat);
      return cp;
    }



    template <int dim>
    std::array<std::pair<double, unsigned int>,2>
    GrainSize<dim>::
    enthalpy_derivative (const typename Interface<dim>::MaterialModelInputs &in) const
    {
      std::array<std::pair<double, unsigned int>,2> derivative;

      if (in.current_cell.state() == IteratorState::valid)
        {
          // get the pressures and temperatures at the vertices of the cell
          const QTrapezoid<dim> quadrature_formula;

          std::vector<double> solution_values(this->get_fe().dofs_per_cell);
          in.current_cell->get_dof_values(this->get_current_linearization_point(),
                                          solution_values.begin(),
                                          solution_values.end());

          // Only create the evaluator the first time we get here
          if (!temperature_evaluator)
            temperature_evaluator
              = std::make_unique<FEPointEvaluation<1,dim>>(this->get_mapping(),
                                                            this->get_fe(),
                                                            update_values,
                                                            this->introspection().component_indices.temperature);
          if (!pressure_evaluator)
            pressure_evaluator
              = std::make_unique<FEPointEvaluation<1,dim>>(this->get_mapping(),
                                                            this->get_fe(),
                                                            update_values,
                                                            this->introspection().component_indices.pressure);


          // Initialize the evaluator for the temperature
          temperature_evaluator->reinit(in.current_cell, quadrature_formula.get_points());
          temperature_evaluator->evaluate(solution_values,
                                          EvaluationFlags::values);

          // Initialize the evaluator for the pressure
          pressure_evaluator->reinit(in.current_cell, quadrature_formula.get_points());
          pressure_evaluator->evaluate(solution_values,
                                       EvaluationFlags::values);

          std::vector<double> temperatures(quadrature_formula.size());
          std::vector<double> pressures(quadrature_formula.size());

          for (unsigned int i=0; i<quadrature_formula.size(); ++i)
            {
              temperatures[i] = temperature_evaluator->get_value(i);
              pressures[i] = pressure_evaluator->get_value(i);
            }

          AssertThrow (material_lookup.size() == 1,
                       ExcMessage("This formalism is only implemented for one material "
                                  "table."));

          // We have to take into account here that the p,T spacing of the table of material properties
          // we use might be on a finer grid than our model. Because of that we compute the enthalpy
          // derivatives by using finite differences that average over the whole temperature and
          // pressure range that is used in this cell. This way we should not miss any phase transformation.
          derivative = material_lookup[0]->enthalpy_derivatives(temperatures,
                                                                pressures,
                                                                max_latent_heat_substeps);
        }

      return derivative;
    }



    template <int dim>
    void
    GrainSize<dim>::
    evaluate(const typename Interface<dim>::MaterialModelInputs &in,
             typename Interface<dim>::MaterialModelOutputs &out) const
    {
      std::vector<double> adiabatic_pressures (in.n_evaluation_points());
      std::vector<unsigned int> phase_indices (in.n_evaluation_points());

      const unsigned int grain_size_index = this->introspection().get_indices_for_fields_of_type(CompositionalFieldDescription::grain_size)[0];

      for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
        {
          // Use the adiabatic pressure instead of the real one, because of oscillations
          adiabatic_pressures[i] = (this->get_adiabatic_conditions().is_initialized())
                                   ?
                                   this->get_adiabatic_conditions().pressure(in.position[i])
                                   :
                                   in.pressure[i];

          const double gravity_norm = this->get_gravity_model().gravity_vector(in.position[i]).norm();

          out.densities[i] = density(in.temperature[i], adiabatic_pressures[i], in.composition[i], in.position[i]);
          out.thermal_conductivities[i] = k_value;
          out.compressibilities[i] = compressibility(in.temperature[i], adiabatic_pressures[i], in.composition[i], in.position[i]);

          // We do not fill the phase function index, because that will be done internally in the get_phase_index() function
          const double depth = this->get_geometry_model().depth(in.position[i]);
          const double rho_g = out.densities[i] * gravity_norm;
          MaterialUtilities::PhaseFunctionInputs<dim> phase_inputs(in.temperature[i], adiabatic_pressures[i], depth, rho_g, numbers::invalid_unsigned_int);
          phase_indices[i] = get_phase_index(phase_inputs);

          if (in.requests_property(MaterialProperties::viscosity))
            {
              double effective_viscosity;
              double disl_viscosity = std::numeric_limits<double>::max();
              Assert(std::isfinite(in.strain_rate[i].norm()),
                     ExcMessage("Invalid strain_rate in the MaterialModelInputs. This is likely because it was "
                                "not filled by the caller."));
              const SymmetricTensor<2,dim> shear_strain_rate = deviator(in.strain_rate[i]);
              const double second_strain_rate_invariant = std::sqrt(std::max(-second_invariant(shear_strain_rate), 0.));

              const double adiabatic_temperature = this->get_adiabatic_conditions().is_initialized()
                                                   ?
                                                   this->get_adiabatic_conditions().temperature(in.position[i])
                                                   :
                                                   in.temperature[i];

              // Make sure grain size is not negative/too small.
              const double limited_grain_size = std::max(minimum_grain_size,in.composition[i][grain_size_index]);
              const double diff_viscosity = diffusion_viscosity(in.temperature[i],
                                                                adiabatic_temperature,
                                                                adiabatic_pressures[i],
                                                                limited_grain_size,
                                                                second_strain_rate_invariant,
                                                                phase_indices[i]);

              if (std::abs(second_strain_rate_invariant) > 1e-30)
                {
                  disl_viscosity = dislocation_viscosity(in.temperature[i], adiabatic_temperature, adiabatic_pressures[i], in.strain_rate[i], phase_indices[i], diff_viscosity);
                  effective_viscosity = disl_viscosity * diff_viscosity / (disl_viscosity + diff_viscosity);
                }
              else
                effective_viscosity = diff_viscosity;

              if (enable_drucker_prager_rheology)
                {
                  // Calculate non-yielding (viscous) stress magnitude.
                  const double non_yielding_stress = 2. * effective_viscosity * second_strain_rate_invariant;

                  // The following handles phases
                  std::vector<unsigned int> n_phases = {n_phase_transitions+1};
                  std::vector<double> phase_function_values(n_phase_transitions, 0.0);

                  for (unsigned int k=0; k<n_phase_transitions; ++k)
                    {
                      phase_inputs.phase_transition_index = k;
                      phase_function_values[k] = phase_function->compute_value(phase_inputs);
                    }

                  // In the grain size material model, viscosity does not depend on composition,
                  // so we set the compositional index for the Drucker-Prager parameters to 0.
                  const Rheology::DruckerPragerParameters drucker_prager_parameters = drucker_prager_plasticity.compute_drucker_prager_parameters(0,
                                                                                      phase_function_values,
                                                                                      n_phases);
                  const double pressure_for_yielding = use_adiabatic_pressure_for_yielding
                                                       ?
                                                       adiabatic_pressures[i]
                                                       :
                                                       std::max(in.pressure[i],0.0);

                  const double yield_stress = drucker_prager_plasticity.compute_yield_stress(drucker_prager_parameters.cohesion,
                                                                                             drucker_prager_parameters.angle_internal_friction,
                                                                                             pressure_for_yielding,
                                                                                             drucker_prager_parameters.max_yield_stress);

                  // Apply plastic yielding:
                  // If the non-yielding stress is greater than the yield stress,
                  // rescale the viscosity back to yield surface
                  if (non_yielding_stress >= yield_stress)
                    {
                      effective_viscosity = drucker_prager_plasticity.compute_viscosity(drucker_prager_parameters.cohesion,
                                                                                        drucker_prager_parameters.angle_internal_friction,
                                                                                        pressure_for_yielding,
                                                                                        second_strain_rate_invariant,
                                                                                        drucker_prager_parameters.max_yield_stress,
                                                                                        effective_viscosity);
                    }

                  PlasticAdditionalOutputs<dim> *plastic_out = out.template get_additional_output<PlasticAdditionalOutputs<dim>>();

                  if (plastic_out != nullptr)
                    {
                      plastic_out->cohesions[i] = drucker_prager_parameters.cohesion;
                      plastic_out->friction_angles[i] = drucker_prager_parameters.angle_internal_friction;
                      plastic_out->yield_stresses[i] = yield_stress;
                      plastic_out->yielding[i] = non_yielding_stress >= yield_stress ? 1 : 0;
                    }
                }

              out.viscosities[i] = std::min(std::max(min_eta,effective_viscosity),max_eta);

              if (DislocationViscosityOutputs<dim> *disl_viscosities_out = out.template get_additional_output<DislocationViscosityOutputs<dim>>())
                {
                  disl_viscosities_out->dislocation_viscosities[i] = std::min(std::max(min_eta,disl_viscosity),1e300);
                  disl_viscosities_out->diffusion_viscosities[i] = std::min(std::max(min_eta,diff_viscosity),1e300);
                }

            }

          // fill seismic velocities outputs if they exist
          if (use_table_properties)
            if (SeismicAdditionalOutputs<dim> *seismic_out = out.template get_additional_output<SeismicAdditionalOutputs<dim>>())
              {
                seismic_out->vp[i] = seismic_Vp(in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
                seismic_out->vs[i] = seismic_Vs(in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
              }
        }

      DislocationViscosityOutputs<dim> *disl_viscosities_out = out.template get_additional_output<DislocationViscosityOutputs<dim>>();
      grain_size_evolution->fill_additional_outputs(in,out,phase_indices,disl_viscosities_out->dislocation_viscosities,out.additional_outputs);

      if (in.requests_property(MaterialProperties::reaction_terms))
        {
          // Create the two lambda functions that are needed to calculate the reaction terms.
          // The functions give access to the dislocation and diffusion viscosity functions of this class.
          const std::function<double(double, double, double, const dealii::SymmetricTensor<2, dim>&, unsigned int, double, double)> dislocation_viscosity_ = [this] (const double temperature,
              const double adiabatic_temperature,
              const double adiabatic_pressure,
              const SymmetricTensor<2,dim> &strain_rate,
              const unsigned int phase_index,
              const double diffusion_viscosity,
              const double viscosity_guess)->double
          {
            return this->dislocation_viscosity(temperature, adiabatic_temperature, adiabatic_pressure, strain_rate, phase_index, diffusion_viscosity, viscosity_guess);
          };

          const std::function<double(double, double, double, double, double, unsigned int)> diffusion_viscosity_ = [this] (const double temperature,
              const double adiabatic_temperature,
              const double adiabatic_pressure,
              const double grain_size,
              const double second_strain_rate_invariant,
              const unsigned int phase_index)->double
          {
            return this->diffusion_viscosity(temperature, adiabatic_temperature, adiabatic_pressure, grain_size, second_strain_rate_invariant, phase_index);
          };

          // Initialize reaction terms.
          for (auto &reaction_term: out.reaction_terms)
            for (auto &reaction_composition: reaction_term)
              reaction_composition = 0.0;

          // Let the grain size evolution model calculate the reaction terms.
          grain_size_evolution->calculate_reaction_terms(in, adiabatic_pressures, phase_indices, dislocation_viscosity_, diffusion_viscosity_, min_eta, max_eta, out);
        }

      /* We separate the calculation of specific heat and thermal expansivity,
       * because they depend on cell-wise averaged values that are only available
       * here
       */
      double average_temperature(0.0);
      double average_density(0.0);
      for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
        {
          average_temperature += in.temperature[i];
          average_density += out.densities[i];
        }
      average_temperature /= in.n_evaluation_points();
      average_density /= in.n_evaluation_points();

      std::array<std::pair<double, unsigned int>,2> dH;

      if (use_table_properties && use_enthalpy)
        dH = enthalpy_derivative(in);

      for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
        {
          if (!use_table_properties)
            {
              out.thermal_expansion_coefficients[i] = thermal_alpha;
              out.specific_heat[i] = reference_specific_heat;
            }
          else if (use_enthalpy)
            {
              if (this->get_adiabatic_conditions().is_initialized()
                  && (in.current_cell.state() == IteratorState::valid)
                  && (dH[0].second > 0)
                  && (dH[1].second > 0))
                {
                  out.thermal_expansion_coefficients[i] = (1 - average_density * dH[1].first) / average_temperature;
                  out.specific_heat[i] = dH[0].first;
                }
              else
                {
                  if (material_lookup.size() == 1)
                    {
                      out.thermal_expansion_coefficients[i] = (1 - out.densities[i] * material_lookup[0]->dHdp(in.temperature[i],adiabatic_pressures[i])) / in.temperature[i];
                      out.specific_heat[i] = material_lookup[0]->dHdT(in.temperature[i],adiabatic_pressures[i]);
                    }
                  else
                    {
                      ExcNotImplemented();
                    }
                }
            }
          else
            {
              out.thermal_expansion_coefficients[i] = thermal_expansion_coefficient(in.temperature[i], adiabatic_pressures[i], in.composition[i], in.position[i]);
              out.specific_heat[i] = specific_heat(in.temperature[i], adiabatic_pressures[i], in.composition[i], in.position[i]);
            }

          out.thermal_expansion_coefficients[i] = std::max(std::min(out.thermal_expansion_coefficients[i],max_thermal_expansivity),min_thermal_expansivity);
          out.specific_heat[i] = std::max(std::min(out.specific_heat[i],max_specific_heat),min_specific_heat);
        }
    }



    template <int dim>
    void
    GrainSize<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Grain size model");
        {
          prm.declare_entry ("Reference density", "3300",
                             Patterns::Double (0.),
                             "The reference density $\\rho_0$. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
          prm.declare_entry ("Reference temperature", "293.",
                             Patterns::Double (0.),
                             "The reference temperature $T_0$. Units: \\si{\\kelvin}.");
          prm.declare_entry ("Viscosity", "5e24",
                             Patterns::Double (0.),
                             "The value of the constant viscosity. "
                             "Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0.),
                             "The value of the thermal conductivity $k$. "
                             "Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");
          prm.declare_entry ("Reference specific heat", "1250.",
                             Patterns::Double (0.),
                             "The value of the specific heat $cp$. "
                             "Units: \\si{\\joule\\per\\kelvin\\per\\kilogram}.");
          prm.declare_entry ("Thermal expansion coefficient", "2e-5",
                             Patterns::Double (0.),
                             "The value of the thermal expansion coefficient $\\alpha$. "
                             "Units: \\si{\\per\\kelvin}.");
          prm.declare_entry ("Reference compressibility", "4e-12",
                             Patterns::Double (0.),
                             "The value of the reference compressibility. "
                             "Units: \\si{\\per\\pascal}.");

          MaterialUtilities::PhaseFunction<dim>::declare_parameters(prm);

          prm.declare_entry ("Dislocation viscosity iteration threshold", "1e-3",
                             Patterns::Double (0.),
                             "We need to perform an iteration inside the computation "
                             "of the dislocation viscosity, because it depends on the "
                             "dislocation strain rate, which depends on the dislocation "
                             "viscosity itself. This number determines the termination "
                             "accuracy, i.e. if the dislocation viscosity changes by less "
                             "than this factor we terminate the iteration.");
          prm.declare_entry ("Dislocation viscosity iteration number", "100",
                             Patterns::Integer(0),
                             "We need to perform an iteration inside the computation "
                             "of the dislocation viscosity, because it depends on the "
                             "dislocation strain rate, which depends on the dislocation "
                             "viscosity itself. This number determines the maximum "
                             "number of iterations that are performed. ");
          prm.declare_entry ("Dislocation creep exponent", "3.5",
                             Patterns::List (Patterns::Double (0.)),
                             "The power-law exponent $n_{dis}$ for dislocation creep. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: none.");
          prm.declare_entry ("Dislocation activation energy", "4.8e5",
                             Patterns::List (Patterns::Double (0.)),
                             "The activation energy for dislocation creep $E_{dis}$. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: \\si{\\joule\\per\\mole}.");
          prm.declare_entry ("Dislocation activation volume", "1.1e-5",
                             Patterns::List (Patterns::Double (0.)),
                             "The activation volume for dislocation creep $V_{dis}$. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: \\si{\\meter\\cubed\\per\\mole}.");
          prm.declare_entry ("Dislocation creep prefactor", "4.5e-15",
                             Patterns::List (Patterns::Double (0.)),
                             "The prefactor for the dislocation creep law $A_{dis}$. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: \\si{\\pascal}$^{-n_{dis}}$\\si{\\per\\second}.");
          prm.declare_entry ("Diffusion creep exponent", "1.",
                             Patterns::List (Patterns::Double (0.)),
                             "The power-law exponent $n_{diff}$ for diffusion creep. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: none.");
          prm.declare_entry ("Diffusion activation energy", "3.35e5",
                             Patterns::List (Patterns::Double (0.)),
                             "The activation energy for diffusion creep $E_{diff}$. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: \\si{\\joule\\per\\mole}.");
          prm.declare_entry ("Diffusion activation volume", "4e-6",
                             Patterns::List (Patterns::Double (0.)),
                             "The activation volume for diffusion creep $V_{diff}$. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: \\si{\\meter\\cubed\\per\\mole}.");
          prm.declare_entry ("Diffusion creep prefactor", "7.4e-15",
                             Patterns::List (Patterns::Double (0.)),
                             "The prefactor for the diffusion creep law $A_{diff}$. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: \\si{\\meter}$^{p_{diff}}$\\si{\\pascal}$^{-n_{diff}}$\\si{\\per\\second}.");
          prm.declare_entry ("Diffusion creep grain size exponent", "3.",
                             Patterns::List (Patterns::Double (0.)),
                             "The diffusion creep grain size exponent $p_{diff}$ that determines the "
                             "dependence of viscosity on grain size. "
                             "List must have one more entry than the Phase transition depths. "
                             "Units: none.");
          prm.declare_entry ("Maximum temperature dependence of viscosity", "100.",
                             Patterns::Double (0.),
                             "The factor by which viscosity at adiabatic temperature and ambient temperature "
                             "are allowed to differ (a value of x means that the viscosity can be x times higher "
                             "or x times lower compared to the value at adiabatic temperature. This parameter "
                             "is introduced to limit local viscosity contrasts, but still allow for a widely "
                             "varying viscosity over the whole mantle range. "
                             "Units: none.");
          prm.declare_entry ("Minimum viscosity", "1e18",
                             Patterns::Double (0.),
                             "The minimum viscosity that is allowed in the whole model domain. "
                             "Units: Pa \\, s.");
          prm.declare_entry ("Maximum viscosity", "1e26",
                             Patterns::Double (0.),
                             "The maximum viscosity that is allowed in the whole model domain. "
                             "Units: Pa \\, s.");
          prm.declare_entry ("Minimum specific heat", "500.",
                             Patterns::Double (0.),
                             "The minimum specific heat that is allowed in the whole model domain. "
                             "Units: J/kg/K.");
          prm.declare_entry ("Maximum specific heat", "6000.",
                             Patterns::Double (0.),
                             "The maximum specific heat that is allowed in the whole model domain. "
                             "Units: J/kg/K.");
          prm.declare_entry ("Minimum thermal expansivity", "1e-5",
                             Patterns::Double (),
                             "The minimum thermal expansivity that is allowed in the whole model domain. "
                             "Units: 1/K.");
          prm.declare_entry ("Maximum thermal expansivity", "1e-3",
                             Patterns::Double (),
                             "The maximum thermal expansivity that is allowed in the whole model domain. "
                             "Units: 1/K.");
          prm.declare_entry ("Maximum latent heat substeps", "1",
                             Patterns::Integer (1),
                             "The maximum number of substeps over the temperature pressure range "
                             "to calculate the averaged enthalpy gradient over a cell.");
          prm.declare_entry ("Minimum grain size", "1e-5",
                             Patterns::Double (0.),
                             "The minimum grain size that is used for the material model. This parameter "
                             "is introduced to limit local viscosity contrasts, but still allows for a widely "
                             "varying viscosity over the whole mantle range. "
                             "Units: \\si{\\meter}.");
          prm.declare_entry ("Lower mantle grain size scaling", "1.0",
                             Patterns::Double (0.),
                             "This option does not exist any more.");
          prm.declare_entry ("Advect logarithm of grain size", "false",
                             Patterns::Bool (),
                             "This option does not exist any more.");
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/material-model/steinberger/",
                             Patterns::DirectoryName (),
                             "The path to the model data. The path may also include the special "
                             "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the 'data/' subdirectory of ASPECT. ");
          prm.declare_entry ("Material file names", "pyr-ringwood88.txt",
                             Patterns::List (Patterns::Anything()),
                             "The file names of the material data. "
                             "List with as many components as active "
                             "compositional fields (material data is assumed to "
                             "be in order with the ordering of the fields). ");
          prm.declare_entry ("Derivatives file names", "",
                             Patterns::List (Patterns::Anything()),
                             "The file names of the enthalpy derivatives data. "
                             "List with as many components as active "
                             "compositional fields (material data is assumed to "
                             "be in order with the ordering of the fields). ");
          prm.declare_entry ("Use table properties", "false",
                             Patterns::Bool(),
                             "This parameter determines whether to use the table properties "
                             "also for density, thermal expansivity and specific heat. "
                             "If false the properties are generated as in the "
                             "simple compressible plugin.");
          prm.declare_entry ("Material file format", "perplex",
                             Patterns::Selection ("perplex|hefesto"),
                             "The material file format to be read in the property "
                             "tables.");
          prm.declare_entry ("Use enthalpy for material properties", "true",
                             Patterns::Bool(),
                             "This parameter determines whether to use the enthalpy to calculate "
                             "the thermal expansivity and specific heat (if true) or use the "
                             "thermal expansivity and specific heat values from "
                             "the material properties table directly (if false).");
          prm.declare_entry ("Bilinear interpolation", "true",
                             Patterns::Bool (),
                             "This parameter determines whether to use bilinear interpolation "
                             "to compute material properties (slower but more accurate).");

          // Drucker Prager plasticity parameters
          prm.declare_entry ("Use Drucker-Prager rheology", "false",
                             Patterns::Bool(),
                             "This parameter determines whether to apply plastic yielding "
                             "according to a Drucker-Prager rheology after computing the viscosity "
                             "from the (grain-size dependent) visous creep flow laws (if true) "
                             "or not (if false).");
          prm.declare_entry ("Use adiabatic pressure for yield stress", "false",
                             Patterns::Bool (),
                             "Whether to use the adiabatic pressure (if true) instead of the full "
                             "(non-negative) pressure (if false) when calculating the yield stress. "
                             "Using the adiabatic pressure (which is analogous to the depth-dependent "
                             "von Mises model) can be useful to avoid the strong non-linearity associated "
                             "with dynamic pressure variations affecting the yield strength, which can "
                             "make the problem ill-posed. However, dynamic pressure can affect the "
                             "localization of the strain rate and the resulting deformation, and neglecting "
                             "it therefore changes the solution.");
          Rheology::DruckerPrager<dim>::declare_parameters(prm);
          ReactionModel::GrainSizeEvolution<dim>::declare_parameters(prm);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    GrainSize<dim>::parse_parameters (ParameterHandler &prm)
    {
      AssertThrow (this->introspection().get_number_of_fields_of_type(CompositionalFieldDescription::grain_size) == 1,
                   ExcMessage("The 'grain size' material model only works if exactly one compositional "
                              "field with type 'grain size' is present. It looks like there are " +
                              std::to_string(this->introspection().get_number_of_fields_of_type(CompositionalFieldDescription::grain_size))
                              + " fields of this type."));

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Grain size model");
        {
          reference_rho              = prm.get_double ("Reference density");
          reference_T                = prm.get_double ("Reference temperature");
          eta                        = prm.get_double ("Viscosity");
          k_value                    = prm.get_double ("Thermal conductivity");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          thermal_alpha              = prm.get_double ("Thermal expansion coefficient");
          reference_compressibility  = prm.get_double ("Reference compressibility");

          // Phase transition parameters
          phase_function = std::make_shared<MaterialUtilities::PhaseFunction<dim>>();
          phase_function->initialize_simulator (this->get_simulator());
          phase_function->parse_parameters (prm);

          std::vector<unsigned int> n_phases_for_each_composition = phase_function->n_phases_for_each_composition();
          n_phase_transitions = n_phases_for_each_composition[0] - 1;

          for (unsigned int i=1; i<n_phase_transitions; ++i)
            AssertThrow(phase_function->get_transition_depth(i-1)<phase_function->get_transition_depth(i),
                        ExcMessage("Error: Phase transition depths have to be sorted in ascending order!"));

          // rheology parameters
          dislocation_viscosity_iteration_threshold = prm.get_double("Dislocation viscosity iteration threshold");
          dislocation_viscosity_iteration_number = prm.get_integer("Dislocation viscosity iteration number");
          dislocation_creep_exponent            = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Dislocation creep exponent")));
          dislocation_activation_energy         = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Dislocation activation energy")));
          dislocation_activation_volume         = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Dislocation activation volume")));
          dislocation_creep_prefactor           = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Dislocation creep prefactor")));
          diffusion_creep_exponent              = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Diffusion creep exponent")));
          diffusion_activation_energy           = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Diffusion activation energy")));
          diffusion_activation_volume           = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Diffusion activation volume")));
          diffusion_creep_prefactor             = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Diffusion creep prefactor")));
          diffusion_creep_grain_size_exponent   = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Diffusion creep grain size exponent")));
          max_temperature_dependence_of_eta     = prm.get_double ("Maximum temperature dependence of viscosity");
          min_eta                               = prm.get_double ("Minimum viscosity");
          max_eta                               = prm.get_double ("Maximum viscosity");
          min_specific_heat                     = prm.get_double ("Minimum specific heat");
          max_specific_heat                     = prm.get_double ("Maximum specific heat");
          min_thermal_expansivity               = prm.get_double ("Minimum thermal expansivity");
          max_thermal_expansivity               = prm.get_double ("Maximum thermal expansivity");
          max_latent_heat_substeps              = prm.get_integer ("Maximum latent heat substeps");
          minimum_grain_size                    = prm.get_double ("Minimum grain size");

          // scale recrystallized grain size, diffusion creep and grain growth prefactor accordingly
          diffusion_creep_prefactor[diffusion_creep_prefactor.size()-1] *= std::pow(1.0,diffusion_creep_grain_size_exponent[diffusion_creep_grain_size_exponent.size()-1]);

          // prefactors never appear without their exponents. perform some calculations here to save time later
          for (unsigned int i=0; i<diffusion_creep_prefactor.size(); ++i)
            diffusion_creep_prefactor[i] = std::pow(diffusion_creep_prefactor[i],-1.0/diffusion_creep_exponent[i]);
          for (unsigned int i=0; i<dislocation_creep_prefactor.size(); ++i)
            dislocation_creep_prefactor[i] = std::pow(dislocation_creep_prefactor[i],-1.0/dislocation_creep_exponent[i]);

          if (dislocation_creep_exponent.size() != dislocation_activation_energy.size() ||
              dislocation_creep_exponent.size() != dislocation_activation_volume.size() ||
              dislocation_creep_exponent.size() != dislocation_creep_prefactor.size() ||
              dislocation_creep_exponent.size() != diffusion_creep_exponent.size() ||
              dislocation_creep_exponent.size() != diffusion_activation_energy.size() ||
              dislocation_creep_exponent.size() != diffusion_activation_volume.size() ||
              dislocation_creep_exponent.size() != diffusion_creep_prefactor.size() ||
              dislocation_creep_exponent.size() != diffusion_creep_grain_size_exponent.size() )
            AssertThrow(false,
                        ExcMessage("Error: The lists of grain size evolution and flow law parameters "
                                   "need to have the same length!"));

          // parameters for reading in tables with material properties
          datadirectory        = prm.get ("Data directory");
          datadirectory = Utilities::expand_ASPECT_SOURCE_DIR(datadirectory);
          material_file_names  = Utilities::split_string_list
                                 (prm.get ("Material file names"));
          derivatives_file_names = Utilities::split_string_list
                                   (prm.get ("Derivatives file names"));
          use_table_properties = prm.get_bool ("Use table properties");
          use_enthalpy = prm.get_bool ("Use enthalpy for material properties");

          // Make sure the grain size field comes after all potential material
          // data fields. Otherwise our material model calculation uses the
          // wrong compositional fields.
          if (use_table_properties && material_file_names.size() > 1)
            {
              AssertThrow(this->introspection().get_indices_for_fields_of_type(CompositionalFieldDescription::grain_size)[0] >= material_file_names.size(),
                          ExcMessage("The compositional fields indicating the major element composition need to be first in the "
                                     "list of compositional fields, but the grain size field seems to have a lower index than the number "
                                     "of provided data files. This is likely inconsistent. Please check the number of provided data "
                                     "files and the order of compositional fields."));
            }

          if (prm.get ("Material file format") == "perplex")
            material_file_format = perplex;
          else if (prm.get ("Material file format") == "hefesto")
            material_file_format = hefesto;
          else
            AssertThrow (false, ExcNotImplemented());

          use_bilinear_interpolation = prm.get_bool ("Bilinear interpolation");

          // Parse plasticity parameters
          enable_drucker_prager_rheology = prm.get_bool ("Use Drucker-Prager rheology");
          use_adiabatic_pressure_for_yielding = prm.get_bool ("Use adiabatic pressure for yield stress");
          drucker_prager_plasticity.initialize_simulator (this->get_simulator());

          std::vector<unsigned int> n_phases = {n_phase_transitions+1};
          drucker_prager_plasticity.parse_parameters(prm, std::make_unique<std::vector<unsigned int>> (n_phases));

          // Parse grain size evolution parameters
          grain_size_evolution = std::make_unique<ReactionModel::GrainSizeEvolution<dim>>();
          grain_size_evolution->initialize_simulator(this->get_simulator());
          grain_size_evolution->initialize_phase_function(phase_function);
          grain_size_evolution->parse_parameters(prm);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();


      // Declare dependencies on solution variables
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;

      this->model_dependence.viscosity = NonlinearDependence::temperature
                                         | NonlinearDependence::pressure
                                         | NonlinearDependence::strain_rate
                                         | NonlinearDependence::compositional_fields;

      this->model_dependence.density = NonlinearDependence::none;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;

      if (use_table_properties)
        {
          this->model_dependence.density |= NonlinearDependence::temperature
                                            | NonlinearDependence::pressure
                                            | NonlinearDependence::compositional_fields;
          this->model_dependence.compressibility = NonlinearDependence::temperature
                                                   | NonlinearDependence::pressure
                                                   | NonlinearDependence::compositional_fields;
          this->model_dependence.specific_heat = NonlinearDependence::temperature
                                                 | NonlinearDependence::pressure
                                                 | NonlinearDependence::compositional_fields;
        }
      else
        {
          if (thermal_alpha != 0)
            this->model_dependence.density |=NonlinearDependence::temperature;
          if (reference_compressibility != 0)
            this->model_dependence.density |=NonlinearDependence::pressure;
        }
    }



    template <int dim>
    void
    GrainSize<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // These properties are useful as output.
      if (out.template get_additional_output<DislocationViscosityOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::DislocationViscosityOutputs<dim>> (n_points));
        }

      // Let the reaction model create additional outputs
      grain_size_evolution->create_additional_named_outputs(out);

      // These properties are only output properties.
      if (use_table_properties && out.template get_additional_output<SeismicAdditionalOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::SeismicAdditionalOutputs<dim>> (n_points));
        }

      if (enable_drucker_prager_rheology && out.template get_additional_output<PlasticAdditionalOutputs<dim>>() == nullptr)
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
    ASPECT_REGISTER_MATERIAL_MODEL(GrainSize,
                                   "grain size",
                                   "A material model that relies on compositional "
                                   "fields that correspond to the average grain sizes of a "
                                   "mineral phase and source terms that determine the grain "
                                   "size evolution in terms of the strain rate, "
                                   "temperature, phase transitions, and the creep regime. "
                                   "This material model only works if a compositional field "
                                   "named 'grain_size' is present. "
                                   "In the diffusion creep regime, the viscosity depends "
                                   "on this grain size field. "
                                   "We use the grain size evolution laws described in Behn "
                                   "et al., 2009. Implications of grain size evolution on the "
                                   "seismic structure of the oceanic upper mantle, "
                                   "Earth Planet. Sci. Letters, 282, 178–189. "
                                   "Other material parameters are either prescribed similar "
                                   "to the 'simple' material model, or read from data files "
                                   "that were generated by the Perplex or Hefesto software. "
                                   "This material model "
                                   "is described in more detail in Dannberg, J., Z. Eilon, "
                                   "U. Faul, R. Gassmoeller, P. Moulik, and R. Myhill (2017), "
                                   "The importance of grain size to mantle dynamics and "
                                   "seismological observations, Geochem. Geophys. Geosyst., "
                                   "18, 3034–3061, doi:10.1002/2017GC006944.")

#define INSTANTIATE(dim) \
  template class DislocationViscosityOutputs<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
