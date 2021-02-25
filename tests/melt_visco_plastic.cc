/*
  Copyright (C) 2015 - 2021 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_melt_visco_plastic_h
#define _aspect_material_model_melt_visco_plastic_h

#include <aspect/material_model/interface.h>
#include <aspect/material_model/visco_plastic.h>
#include <aspect/material_model/equation_of_state/multicomponent_incompressible.h>

#include <aspect/simulator.h>
#include <aspect/simulator_access.h>
#include <aspect/postprocess/melt_statistics.h>
#include <aspect/melt.h>
#include <aspect/utilities.h>
#include <aspect/newton.h>

#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model that implements a simple formulation of the
     * material parameters required for the modelling of melt transport,
     * including a source term for the porosity according to the melting
     * model for dry peridotite of Katz, 2003. This also includes a
     * computation of the latent heat of melting (if the latent heat
     * heating model is active).
     *
     * Most of the material properties are constant, except for the shear,
     * compaction and melt viscosities and the permeability, which depend on
     * the porosity; and the solid and melt densities, which depend on
     * temperature and pressure.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class MeltViscoPlastic : public MaterialModel::MeltFractionModel<dim>,
      public MaterialModel::MeltInterface<dim>,
      public ::aspect::SimulatorAccess<dim>
    {
      public:

        /**
         * Return whether the model is compressible or not.
         * In this case, Keller et al, 2003 assume the fluid and
         * solid phases are each incompressible: "compressibility
         * in the model is accounted for by changes in melt fraction."
         * Thus, this function will always return false.
         */
        virtual bool is_compressible () const override;

        /**
         * @name Reference quantities
         * @{
         */
        virtual double reference_viscosity () const override;

        virtual double reference_darcy_coefficient () const override;

        virtual void evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                              typename Interface<dim>::MaterialModelOutputs &out) const override;

        virtual void melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                                     std::vector<double> &melt_fractions) const;

        /**
         * @}
         */

        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */
        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);
        /**
         * @}
         */

        virtual
        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const override;

      private:

        double min_strain_rate;
        double ref_strain_rate;
        double min_viscosity;
        double max_viscosity;

        double ref_viscosity;

        /**
         * Enumeration for selecting which viscosity averaging scheme to use.
         */

        std::vector<double> linear_viscosities;

        std::vector<double> angles_internal_friction;
        std::vector<double> cohesions;

        double maximum_yield_stress;

        std::vector<double> strength_reductions;

        double melt_density_change;
        double xi_0;
        double eta_f;
        double reference_permeability;
        double alpha_phi;
        double freezing_rate;
        double melting_time_scale;

        /**
         * Parameters for anhydrous melting of peridotite after Katz, 2003
         */

        // for the solidus temperature
        double A1;   // °C
        double A2; // °C/Pa
        double A3; // °C/(Pa^2)

        // for the lherzolite liquidus temperature
        double B1;   // °C
        double B2;   // °C/Pa
        double B3; // °C/(Pa^2)

        // for the liquidus temperature
        double C1;   // °C
        double C2;  // °C/Pa
        double C3; // °C/(Pa^2)

        // for the reaction coefficient of pyroxene
        double r1;     // cpx/melt
        double r2;     // cpx/melt/GPa
        double M_cpx;  // mass fraction of pyroxene

        // melt fraction exponent
        double beta;

        /**
         * Percentage of material that is molten for a given temperature and
         * pressure (assuming equilibrium conditions). Melting model after Katz,
         * 2003, for dry peridotite.
         */
        virtual
        double
        melt_fraction (const double temperature,
                       const double pressure) const;

        /**
         * A function that fills the plastic additional output in the
         * MaterialModelOutputs object that is handed over, if it exists.
         * Does nothing otherwise.
         */
        // void fill_plastic_outputs (const unsigned int point_index,
        //                            const std::vector<double> &volume_fractions,
        //                            const bool plastic_yielding,
        //                            const MaterialModel::MaterialModelInputs<dim> &in,
        //                            MaterialModel::MaterialModelOutputs<dim> &out) const;

        /**
         * Whether to use a function defined as an input parameter for the melting rate
         */
        bool use_melting_rate_function;

        /**
         * A function object representing the melting rate.
         */
        Functions::ParsedFunction<dim> function;


        /**
         * Pointer to the object used to compute the rheological properties.
         * In this case, the rheology in question is visco(elasto)plastic. The
         * object contains functions for parameter declaration and parsing,
         * and further functions that calculate viscosity and viscosity
         * derivatives. It also contains functions that create and fill
         * additional material model outputs, specifically plastic outputs.
         * The rheology itself is a composite rheology, and so the object
         * contains further objects and/or pointers to objects that provide
         * functions and parameters for all subordinate rheologies.
         */
        std::unique_ptr<Rheology::ViscoPlastic<dim>> rheology;

        std::vector<double> thermal_diffusivities;

        /**
         * Whether to use user-defined thermal conductivities instead of thermal diffusivities.
         */
        bool define_conductivities;

        std::vector<double> thermal_conductivities;

        /**
         * Object for computing the equation of state.
         */
        EquationOfState::MulticomponentIncompressible<dim> equation_of_state;

        /**
         * Object that handles phase transitions.
         */
        MaterialUtilities::PhaseFunction<dim> phase_function;

    };

  }
}

#endif

namespace aspect
{
  namespace MaterialModel
  {

    using namespace dealii;

    template <int dim>
    double
    MeltViscoPlastic<dim>::
    reference_viscosity () const
    {
      return rheology->ref_visc;
    }

    template <int dim>
    double
    MeltViscoPlastic<dim>::
    reference_darcy_coefficient () const
    {
      // 0.01 = 1% melt
      return reference_permeability * std::pow(0.01,3.0) / eta_f;
    }

    template <int dim>
    bool
    MeltViscoPlastic<dim>::
    is_compressible () const
    {
      return equation_of_state.is_compressible();
      // return false;
    }

    template <int dim>
    double
    MeltViscoPlastic<dim>::
    melt_fraction (const double temperature,
                   const double pressure) const
    {
      // anhydrous melting of peridotite after Katz, 2003
      const double T_solidus  = A1 + 273.15
                                + A2 * pressure
                                + A3 * pressure * pressure;
      const double T_lherz_liquidus = B1 + 273.15
                                      + B2 * pressure
                                      + B3 * pressure * pressure;
      const double T_liquidus = C1 + 273.15
                                + C2 * pressure
                                + C3 * pressure * pressure;

      // melt fraction for peridotite with clinopyroxene
      double peridotite_melt_fraction;
      if (temperature < T_solidus || pressure > 1.3e10)
        peridotite_melt_fraction = 0.0;
      else if (temperature > T_lherz_liquidus)
        peridotite_melt_fraction = 1.0;
      else
        peridotite_melt_fraction = std::pow((temperature - T_solidus) / (T_lherz_liquidus - T_solidus),beta);

      // melt fraction after melting of all clinopyroxene
      const double R_cpx = r1 + r2 * std::max(0.0, pressure);
      const double F_max = M_cpx / R_cpx;

      if (peridotite_melt_fraction > F_max && temperature < T_liquidus)
        {
          const double T_max = std::pow(F_max,1/beta) * (T_lherz_liquidus - T_solidus) + T_solidus;
          peridotite_melt_fraction = F_max + (1 - F_max) * pow((temperature - T_max) / (T_liquidus - T_max),beta);
        }
      return peridotite_melt_fraction;
    }

    template <int dim>
    void
    MeltViscoPlastic<dim>::
    melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                    std::vector<double> &melt_fractions) const
    {
      for (unsigned int q=0; q<in.temperature.size(); ++q)
        melt_fractions[q] = melt_fraction(in.temperature[q],
                                          std::max(0.0, in.pressure[q]));
      return;
    }

    template <int dim>
    void
    MeltViscoPlastic<dim>::
    evaluate(const typename Interface<dim>::MaterialModelInputs &in,
             typename Interface<dim>::MaterialModelOutputs &out) const
    {
      // Store which components do not represent volumetric compositions (e.g. strain components).
      const ComponentMask volumetric_compositions = rheology->get_volumetric_composition_mask();

      EquationOfStateOutputs<dim> eos_outputs (this->n_compositional_fields()+1);
      EquationOfStateOutputs<dim> eos_outputs_all_phases (this->n_compositional_fields()+1+phase_function.n_phase_transitions());

      std::vector<double> average_elastic_shear_moduli (in.n_evaluation_points());

      // Store value of phase function for each phase and composition
      // While the number of phases is fixed, the value of the phase function is updated for every point
      std::vector<double> phase_function_values(phase_function.n_phase_transitions(), 0.0);

      // 1) Initial viscosities and other material properties
      for (unsigned int i=0; i<in.position.size(); ++i)
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
          for (unsigned int j=0; j < phase_function.n_phase_transitions(); j++)
            {
              phase_inputs.phase_index = j;
              phase_function_values[j] = phase_function.compute_value(phase_inputs);
            }

          // Average by value of gamma function to get value of compositions
          phase_average_equation_of_state_outputs(eos_outputs_all_phases,
                                                  phase_function_values,
                                                  phase_function.n_phase_transitions_for_each_composition(),
                                                  eos_outputs);

          const std::vector<double> volume_fractions = MaterialUtilities::compute_composition_fractions(in.composition[0], rheology->get_volumetric_composition_mask());
          out.viscosities[i] = MaterialUtilities::average_value(volume_fractions, linear_viscosities, rheology->viscosity_averaging);

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
        }

      // Store the intrinsic (linear) viscosities for computing the compaction viscosities later on
      // (Keller et al. eq. 51).
      const std::vector<double> intrinsic_viscosities = out.viscosities;

      // 2) Retrieve fluid pressure and volumetric strain rate
      std::vector<double> fluid_pressures(in.position.size());
      std::vector<double> volumetric_strain_rates(in.position.size());
      std::vector<double> volumetric_yield_strength(in.position.size());

      ReactionRateOutputs<dim> *reaction_rate_out = out.template get_additional_output<ReactionRateOutputs<dim> >();

      if (this->include_melt_transport() )
        {
          if (in.current_cell.state() == IteratorState::valid)
            {
              std::vector<Point<dim> > quadrature_positions(in.position.size());

              for (unsigned int i=0; i < in.position.size(); ++i)
                quadrature_positions[i] = this->get_mapping().transform_real_to_unit_cell(in.current_cell, in.position[i]);

              // FEValues requires a quadrature and we provide the default quadrature
              // as we only need to evaluate the solution and gradients.
              FEValues<dim> fe_values (this->get_mapping(),
                                       this->get_fe(),
                                       Quadrature<dim>(quadrature_positions),
                                       update_values | update_gradients);

              fe_values.reinit (in.current_cell);

              // get fluid pressure from the current solution
              const FEValuesExtractors::Scalar extractor_pressure = this->introspection().variable("fluid pressure").extractor_scalar();
              fe_values[extractor_pressure].get_function_values (this->get_solution(),
                                                                 fluid_pressures);

              // get volumetric strain rate
              // see Keller et al. eq. 11.
              fe_values[this->introspection().extractors.velocities].get_function_divergences (this->get_solution(),
                  volumetric_strain_rates);
            }

          const double time_scale = this->convert_output_to_years() ? year_in_seconds : 1.0;

          // 3) Get porosity, melt density and update melt reaction terms
          for (unsigned int i=0; i<in.position.size(); ++i)
            {
              // get peridotite and porosity field indices
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
              const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");

              const double old_porosity = in.composition[i][porosity_idx];
              const double maximum_melt_fraction = in.composition[i][peridotite_idx];

              // calculate the melting rate as difference between the equilibrium melt fraction
              // and the solution of the previous time step
              double porosity_change = 0.0;

              // batch melting
              porosity_change = melt_fraction(in.temperature[i], this->get_adiabatic_conditions().pressure(in.position[i]))
                                - std::max(maximum_melt_fraction, 0.0);
              porosity_change = std::max(porosity_change, 0.0);

              // do not allow negative porosity
              porosity_change = std::max(porosity_change, -old_porosity);

              // because depletion is a volume-based, and not a mass-based property that is advected,
              // additional scaling factors on the right hand side apply
              for (unsigned int c=0; c<in.composition[i].size(); ++c)
                {
                  if (use_melting_rate_function == false)
                    {
                      // fill reaction rate outputs
                      if (reaction_rate_out != nullptr)
                        {
                          if (c == peridotite_idx && this->get_timestep_number() > 0)
                            reaction_rate_out->reaction_rates[i][c] = porosity_change / melting_time_scale
                                                                      * (1 - maximum_melt_fraction) / (1 - old_porosity);
                          else if (c == porosity_idx && this->get_timestep_number() > 0)
                            reaction_rate_out->reaction_rates[i][c] = porosity_change / melting_time_scale;
                          else
                            reaction_rate_out->reaction_rates[i][c] = 0.0;
                        }

                      out.reaction_terms[i][c] = 0.0;
                    }
                  else
                    {
                      // Set reaction rates to 0
                      if (reaction_rate_out != nullptr)
                        reaction_rate_out->reaction_rates[i][c] = 0.;
                      //
                      if (c == porosity_idx)
                        out.reaction_terms[i][c] = function.value(in.position[i]) / time_scale * out.densities[i];
                      else if (c == peridotite_idx)
                        out.reaction_terms[i][c] = 0.0;
                    }
                }

              // 4) Reduce shear viscosity due to melt presence
              const double porosity = std::min(1.0, std::max(in.composition[i][porosity_idx],0.0));
              out.viscosities[i] *= exp(- alpha_phi * porosity);
            }
        }

      if (in.strain_rate.size() )
        {
          // 5) Compute plastic weakening of the viscosity
          for (unsigned int i=0; i<in.position.size(); ++i)
            {
              // Compute volume fractions
              const std::vector<double> volume_fractions = MaterialUtilities::compute_composition_fractions(in.composition[i]);

              // // Compute the effective viscosity if requested and retrieve whether the material is plastically yielding
              bool plastic_yielding = false;

              /***************
               * Ideally, the rest of this block can eventually be replaced with the contents of this
               * long comment block (or something very similar). Most of this functionality is included
               * in the `rheology->calculate_isostrain_viscosities` method.

                            if (in.requests_property(MaterialProperties::viscosity))
                              {
                                // Currently, the viscosities for each of the compositional fields are calculated assuming
                                // isostrain amongst all compositions, allowing calculation of the viscosity ratio.
                                // TODO: This is only consistent with viscosity averaging if the arithmetic averaging
                                // scheme is chosen. It would be useful to have a function to calculate isostress viscosities.
                                const IsostrainViscosities isostrain_viscosities =
                                  rheology->calculate_isostrain_viscosities(in, i, volume_fractions, phase_function_values,
                                                                            phase_function.n_phase_transitions_for_each_composition());

                                // The isostrain condition implies that the viscosity averaging should be arithmetic (see above).
                                // We have given the user freedom to apply alternative bounds, because in diffusion-dominated
                                // creep (where n_diff=1) viscosities are stress and strain-rate independent, so the calculation
                                // of compositional field viscosities is consistent with any averaging scheme.
                                out.viscosities[i] = MaterialUtilities::average_value(volume_fractions, isostrain_viscosities.composition_viscosities,
                                                                                      rheology->viscosity_averaging);

                                // Decide based on the maximum composition if material is yielding.
                                // This avoids for example division by zero for harmonic averaging (as plastic_yielding
                                // holds values that are either 0 or 1), but might not be consistent with the viscosity
                                // averaging chosen.
                                std::vector<double>::const_iterator max_composition = std::max_element(volume_fractions.begin(),volume_fractions.end());
                                plastic_yielding = isostrain_viscosities.composition_yielding[std::distance(volume_fractions.begin(),max_composition)];

                                // Compute viscosity derivatives if they are requested
                                if (MaterialModel::MaterialModelDerivatives<dim> *derivatives =
                                      out.template get_additional_output<MaterialModel::MaterialModelDerivatives<dim> >())
                                  rheology->compute_viscosity_derivatives(i, volume_fractions, isostrain_viscosities.composition_viscosities, in, out,
                                                                          phase_function_values, phase_function.n_phase_transitions_for_each_composition());
                              }

                            // Now compute changes in the compositional fields (i.e. the accumulated strain).
                            for (unsigned int c=0; c<in.composition[i].size(); ++c)
                              out.reaction_terms[i][c] = 0.0;

                            // Calculate changes in strain invariants and update the reaction terms
                            rheology->strain_rheology.fill_reaction_outputs(in, i, rheology->min_strain_rate, plastic_yielding, out);

                            // Fill plastic outputs if they exist.
                            rheology->fill_plastic_outputs(i,volume_fractions,plastic_yielding,in,out);

                            if (rheology->use_elasticity)
                              {
                                // Compute average elastic shear modulus
                                average_elastic_shear_moduli[i] = MaterialUtilities::average_value(volume_fractions,
                                                                                                  rheology->elastic_rheology.get_elastic_shear_moduli(),
                                                                                                  rheology->viscosity_averaging);

                                // Fill the material properties that are part of the elastic additional outputs
                                if (ElasticAdditionalOutputs<dim> *elastic_out = out.template get_additional_output<ElasticAdditionalOutputs<dim> >())
                                  {
                                    elastic_out->elastic_shear_moduli[i] = average_elastic_shear_moduli[i];
                                  }
                              }
              ***************/

              // 4) Compute plastic weakening of visco(elastic) viscosity
              double porosity = 0.0;

              if (this->include_melt_transport() )
                porosity = std::min(1.0, std::max(in.composition[i][this->introspection().compositional_index_for_name("porosity")],0.0));

              // calculate deviatoric strain rate (Keller et al. eq. 13)
              const double edot_ii = ( (this->get_timestep_number() == 0 && in.strain_rate[i].norm() <= std::numeric_limits<double>::min())
                                       ?
                                       ref_strain_rate
                                       :
                                       std::max(std::sqrt(std::fabs(second_invariant(deviator(in.strain_rate[i])))),
                                                min_strain_rate) );

              // compute viscous stress
              const double viscous_stress = 2. * out.viscosities[i] * edot_ii * (1.0 - porosity);

              // In case porosity lies above the melt transport threshold
              // P_effective = P_bulk - P_f = (1-porosity) * P_s + porosity * P_f - P_f = (1-porosity) * (P_s - P_f)
              // otherwise,
              // P_effective = P_bulk, which equals P_solid (which is given by in.pressure[i])
              const double effective_pressure = ((this->include_melt_transport() && this->get_melt_handler().is_melt_cell(in.current_cell))
                                                 ?
                                                 (1. - porosity) * (in.pressure[i] - fluid_pressures[i])
                                                 :
                                                 in.pressure[i]);

              double yield_strength = 0.0;
              double tensile_strength = 0.0;

              for (unsigned int c=0; c< volume_fractions.size(); ++c)
                {
                  const double tensile_strength_c = cohesions[c]/strength_reductions[c];

                  // Convert friction angle from degrees to radians
                  double phi = angles_internal_friction[c] * numbers::PI/180.0;
                  const double transition_pressure = (cohesions[c] * std::cos(phi) - tensile_strength_c) / (1.0 -  sin(phi));

                  double yield_strength_c = 0.0;
                  // In case we're not using the Keller et al. formulation,
                  // or the effective pressure is bigger than the transition pressure, use
                  // the normal yield strength formulation
                  if (effective_pressure > transition_pressure || !this->include_melt_transport())
                    yield_strength_c = ( (dim==3)
                                         ?
                                         ( 6.0 * cohesions[c] * std::cos(phi) + 2.0 * effective_pressure * std::sin(phi) )
                                         / ( std::sqrt(3.0) * (3.0 + std::sin(phi) ) )
                                         :
                                         cohesions[c] * std::cos(phi) + effective_pressure * std::sin(phi) );
                  else
                    // Note typo in Keller et al. paper eq. (37) (minus sign)
                    yield_strength_c = tensile_strength_c + effective_pressure;

                  // TODO add different averagings?
                  yield_strength += volume_fractions[c]*yield_strength_c;
                  tensile_strength += volume_fractions[c]*tensile_strength_c;
                }

              // If the viscous stress is greater than the yield strength, rescale the viscosity back to yield surface
              // and reaction term for plastic finite strain
              if (viscous_stress >= yield_strength)
                {
                  plastic_yielding = true;
                  out.viscosities[i] = yield_strength / (2.0 * edot_ii);
                }

              // Calculate average internal friction angle and cohesion values
              const double cohesion = MaterialUtilities::average_value(volume_fractions, cohesions, rheology->viscosity_averaging);
              const double angle_internal_friction = MaterialUtilities::average_value(volume_fractions, angles_internal_friction, rheology->viscosity_averaging);

              PlasticAdditionalOutputs<dim> *plastic_out = out.template get_additional_output<PlasticAdditionalOutputs<dim> >();
              if (plastic_out != nullptr)
                {
                  plastic_out->cohesions[i] = cohesion;
                  plastic_out->friction_angles[i] = angle_internal_friction;
                  plastic_out->yielding[i] = plastic_yielding ? 1 : 0;
                }

              // Limit the viscosity with specified minimum and maximum bounds
              out.viscosities[i] = std::min(std::max(out.viscosities[i], min_viscosity), max_viscosity);

              // Compute the volumetric yield strength (Keller et al. eq (38))
              volumetric_yield_strength[i] = viscous_stress - tensile_strength;
            }
        }

      // fill melt outputs if they exist
      MeltOutputs<dim> *melt_out = out.template get_additional_output<MeltOutputs<dim> >();

      if (melt_out != NULL)
        {
          for (unsigned int i=0; i<in.position.size(); ++i)
            {
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
              double porosity = std::min(1.0, std::max(in.composition[i][porosity_idx],0.0));
              melt_out->fluid_viscosities[i] = eta_f;
              melt_out->permeabilities[i] = reference_permeability * std::pow(porosity,3) * std::pow(1.0-porosity,2);

              melt_out->fluid_densities[i] = out.densities[i] + melt_density_change;
              melt_out->fluid_density_gradients[i] = 0.0;

              const double compaction_pressure = (1.0 - porosity) * (in.pressure[i] - fluid_pressures[i]);

              const double phi_0 = 0.05;
              porosity = std::max(std::min(porosity,0.995),1.e-8);
              // compaction viscosities (Keller et al. eq (51)
              melt_out->compaction_viscosities[i] = intrinsic_viscosities[i] * phi_0 / porosity;

              // visco(elastic) compaction viscosity
              // Keller et al. eq (36)
              // TODO include elastic part
              melt_out->compaction_viscosities[i] *= (1. - porosity);

              // TODO compaction stress evolution parameter
              // Keller et al. eq. (41) and (44)

              // effective compaction viscosity (Keller et al. eq (43) )
              // NB: I've added a minus sign as according to eq 43
              if (in.strain_rate.size() && compaction_pressure < volumetric_yield_strength[i])
                {
                  // the volumetric strain rate might be negative, but will always have the same sign as the volumetric yield strength
                  if (volumetric_strain_rates[i] >= 0)
                    volumetric_strain_rates[i] = std::max(volumetric_strain_rates[i], min_strain_rate);
                  else
                    volumetric_strain_rates[i] = std::min(volumetric_strain_rates[i], -min_strain_rate);
                  melt_out->compaction_viscosities[i] = - volumetric_yield_strength[i] / volumetric_strain_rates[i];
                }

              // Limit the viscosity with specified minimum and maximum bounds
              melt_out->compaction_viscosities[i] = std::min(std::max(melt_out->compaction_viscosities[i], min_viscosity), max_viscosity);
            }
        }
    }


    template <int dim>
    void
    MeltViscoPlastic<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt visco plastic");
        {
          MaterialUtilities::PhaseFunction<dim>::declare_parameters(prm);

          EquationOfState::MulticomponentIncompressible<dim>::declare_parameters (prm);

          Rheology::ViscoPlastic<dim>::declare_parameters(prm);

          // Equation of state parameters
          prm.declare_entry ("Thermal diffusivities", "0.8e-6",
                             Patterns::List(Patterns::Double (0.)),
                             "List of thermal diffusivities, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value.  "
                             "Units: \\si{\\meter\\squared\\per\\second}.");
          prm.declare_entry ("Define thermal conductivities","false",
                             Patterns::Bool (),
                             "Whether to directly define thermal conductivities for each compositional field "
                             "instead of calculating the values through the specified thermal diffusivities, "
                             "densities, and heat capacities. ");

          ////////////////////////

          prm.declare_entry ("Reference temperature", "293", Patterns::Double(0),
                             "For calculating density by thermal expansivity. Units: $K$");
          prm.declare_entry ("Thermal conductivities", "4.7",
                             Patterns::List(Patterns::Double(0)),
                             "List of thermal conductivities for background mantle and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "If only one value is given, then all use the same value. Units: $W/m/K$ ");

          prm.declare_entry ("Minimum strain rate", "1.0e-20", Patterns::Double(0),
                             "Stabilizes strain dependent viscosity. Units: $1 / s$");
          prm.declare_entry ("Reference strain rate","1.0e-15",Patterns::Double(0),
                             "Reference strain rate for first time step. Units: $1 / s$");
          prm.declare_entry ("Minimum viscosity", "1e17", Patterns::Double(0),
                             "Lower cutoff for effective viscosity. Units: $Pa \\, s$");
          prm.declare_entry ("Maximum viscosity", "1e28", Patterns::Double(0),
                             "Upper cutoff for effective viscosity. Units: $Pa \\, s$");

          prm.declare_entry ("Linear viscosities", "1.e22",
                             Patterns::List(Patterns::Double(0)),
                             "List of linear viscosities for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "The values can be used instead of the viscosities derived from the "
                             "base material model. Units: Pa s.");

          prm.declare_entry ("Angles of internal friction", "0",
                             Patterns::List(Patterns::Double(0)),
                             "List of angles of internal friction, $\\phi$, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "For a value of zero, in 2D the von Mises criterion is retrieved. "
                             "Angles higher than 30 degrees are harder to solve numerically. Units: degrees.");
          prm.declare_entry ("Cohesions", "1e20",
                             Patterns::List(Patterns::Double(0)),
                             "List of cohesions, $C$, for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "The extremely large default cohesion value (1e20 Pa) prevents the viscous stress from "
                             "exceeding the yield stress. Units: $Pa$.");

          prm.declare_entry ("Maximum yield stress", "1e12", Patterns::Double(0),
                             "Limits the maximum value of the yield stress determined by the "
                             "drucker-prager plasticity parameters. Default value is chosen so this "
                             "is not automatically used. Values of 100e6--1000e6 $Pa$ have been used "
                             "in previous models. Units: $Pa$");

          prm.declare_entry ("Host rock strength reductions", "4",
                             Patterns::List(Patterns::Double(0)),
                             "List of reduction factors of strength of the host rock under tensile stress, $R$, "
                             "for background material and compositional fields, "
                             "for a total of N+1 values, where N is the number of compositional fields. "
                             "Units: none.");

          prm.declare_entry ("Melt density change", "-500",
                             Patterns::Double (),
                             "Difference between solid density $\\rho_{s}$ and melt/fluid$\\rho_{f}$. "
                             "Units: $kg/m^3$.");
          prm.declare_entry ("Reference bulk viscosity", "1e22",
                             Patterns::Double (0),
                             "The value of the constant bulk viscosity $\\xi_0$ of the solid matrix. "
                             "This viscosity may be modified by both temperature and porosity "
                             "dependencies. Units: $Pa s$.");
          prm.declare_entry ("Reference melt viscosity", "10",
                             Patterns::Double (0),
                             "The value of the constant melt viscosity $\\eta_f$. Units: $Pa s$.");
          prm.declare_entry ("Exponential melt weakening factor", "27",
                             Patterns::Double (0),
                             "The porosity dependence of the viscosity. Units: dimensionless.");
          prm.declare_entry ("Reference permeability", "1e-8",
                             Patterns::Double(),
                             "Reference permeability of the solid host rock."
                             "Units: $m^2$.");
          prm.declare_entry ("Melting time scale for operator splitting", "1e3",
                             Patterns::Double (0),
                             "Because the operator splitting scheme is used, the porosity field can not "
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
                             "computed. "
                             "Units: yr or s, depending on the ``Use years in output instead of seconds'' parameter.");

          prm.declare_entry ("A1", "1085.7",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the solidus "
                             "of peridotite. "
                             "Units: ${}^\\circ C$.");
          prm.declare_entry ("A2", "1.329e-7",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of peridotite. "
                             "Units: ${}^\\circ C/Pa$.");
          prm.declare_entry ("A3", "-5.1e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of peridotite. "
                             "Units: ${}^\\circ C/(Pa^2)$.");
          prm.declare_entry ("B1", "1475.0",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the lherzolite "
                             "liquidus used for calculating the fraction "
                             "of peridotite-derived melt. "
                             "Units: ${}^\\circ C$.");
          prm.declare_entry ("B2", "8.0e-8",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the  lherzolite liquidus used for "
                             "calculating the fraction of peridotite-"
                             "derived melt. "
                             "Units: ${}^\\circ C/Pa$.");
          prm.declare_entry ("B3", "-3.2e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the  lherzolite liquidus used for "
                             "calculating the fraction of peridotite-"
                             "derived melt. "
                             "Units: ${}^\\circ C/(Pa^2)$.");
          prm.declare_entry ("C1", "1780.0",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the liquidus "
                             "of peridotite. "
                             "Units: ${}^\\circ C$.");
          prm.declare_entry ("C2", "4.50e-8",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the liquidus of peridotite. "
                             "Units: ${}^\\circ C/Pa$.");
          prm.declare_entry ("C3", "-2.0e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the liquidus of peridotite. "
                             "Units: ${}^\\circ C/(Pa^2)$.");
          prm.declare_entry ("r1", "0.5",
                             Patterns::Double (),
                             "Constant in the linear function that "
                             "approximates the clinopyroxene reaction "
                             "coefficient. "
                             "Units: non-dimensional.");
          prm.declare_entry ("r2", "8e-11",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the linear function that approximates "
                             "the clinopyroxene reaction coefficient. "
                             "Units: $1/Pa$.");
          prm.declare_entry ("beta", "1.5",
                             Patterns::Double (),
                             "Exponent of the melting temperature in "
                             "the melt fraction calculation. "
                             "Units: non-dimensional.");
          prm.declare_entry ("Mass fraction cpx", "0.15",
                             Patterns::Double (),
                             "Mass fraction of clinopyroxene in the "
                             "peridotite to be molten. "
                             "Units: non-dimensional.");

          prm.declare_entry ("Use melting rate function", "false",
                             Patterns::Bool (),
                             "Prescribe the melting rate through a function defined as an input parameter. "
                             "Units: None");

          prm.enter_subsection("Melting rate function");
          {
            Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
          }
          prm.leave_subsection();

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    MeltViscoPlastic<dim>::parse_parameters (ParameterHandler &prm)
    {
      //increment by one for background:
      const unsigned int n_fields = this->n_compositional_fields() + 1;

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt visco plastic");
        {
          // Phase transition parameters
          phase_function.initialize_simulator (this->get_simulator());
          phase_function.parse_parameters (prm);

          std::vector<unsigned int> n_phase_transitions_for_each_composition
          (phase_function.n_phase_transitions_for_each_composition());

          // We require one more entry for density, etc as there are phase transitions
          // (for the low-pressure phase before any transition).
          for (unsigned int &n : n_phase_transitions_for_each_composition)
            n += 1;

          // Equation of state parameters
          equation_of_state.initialize_simulator (this->get_simulator());
          equation_of_state.parse_parameters (prm,
                                              std::make_shared<std::vector<unsigned int>>(n_phase_transitions_for_each_composition));


          thermal_diffusivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal diffusivities"))),
                                                                          n_fields,
                                                                          "Thermal diffusivities");

          define_conductivities = prm.get_bool ("Define thermal conductivities");

          thermal_conductivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal conductivities"))),
                                                                           n_fields,
                                                                           "Thermal conductivities");

          rheology = std_cxx14::make_unique<Rheology::ViscoPlastic<dim>>();
          rheology->initialize_simulator (this->get_simulator());
          rheology->parse_parameters(prm, std::make_shared<std::vector<unsigned int>>(n_phase_transitions_for_each_composition));

          ////////////////////

          min_strain_rate = prm.get_double("Minimum strain rate");
          ref_strain_rate = prm.get_double("Reference strain rate");
          min_viscosity = prm.get_double ("Minimum viscosity");
          max_viscosity = prm.get_double ("Maximum viscosity");
          ref_viscosity = prm.get_double ("Reference viscosity");

          linear_viscosities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Linear viscosities"))),
                                                                       n_fields,
                                                                       "Linear viscosities");

          angles_internal_friction = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Angles of internal friction"))),
                                                                             n_fields,
                                                                             "Angles of internal friction");
          cohesions = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Cohesions"))),
                                                              n_fields,
                                                              "Cohesions");

          maximum_yield_stress = prm.get_double("Maximum yield stress");

          strength_reductions = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Host rock strength reductions"))),
                                                                        n_fields,
                                                                        "Host rock strength reductions");

          melt_density_change        = prm.get_double ("Melt density change");
          xi_0                       = prm.get_double ("Reference bulk viscosity");
          eta_f                      = prm.get_double ("Reference melt viscosity");
          reference_permeability     = prm.get_double ("Reference permeability");
          alpha_phi                  = prm.get_double ("Exponential melt weakening factor");
          melting_time_scale         = prm.get_double ("Melting time scale for operator splitting");

          A1              = prm.get_double ("A1");
          A2              = prm.get_double ("A2");
          A3              = prm.get_double ("A3");
          B1              = prm.get_double ("B1");
          B2              = prm.get_double ("B2");
          B3              = prm.get_double ("B3");
          C1              = prm.get_double ("C1");
          C2              = prm.get_double ("C2");
          C3              = prm.get_double ("C3");
          r1              = prm.get_double ("r1");
          r2              = prm.get_double ("r2");
          beta            = prm.get_double ("beta");
          M_cpx           = prm.get_double ("Mass fraction cpx");

          use_melting_rate_function  = prm.get_bool ("Use melting rate function");

          if (use_melting_rate_function)
            {
              prm.enter_subsection("Melting rate function");
              try
                {
                  function.parse_parameters (prm);
                }
              catch (...)
                {
                  std::cerr << "ERROR: FunctionParser failed to parse\n"
                            << "\t'Material model.Melting rate function.Function'\n"
                            << "with expression\n"
                            << "\t'" << prm.get("Function expression") << "'"
                            << "More information about the cause of the parse error \n"
                            << "is shown below.\n";
                  throw;
                }
              prm.leave_subsection();
            }

          AssertThrow(this->introspection().compositional_name_exists("peridotite"),
                      ExcMessage("Material model Melt visco plastic only works if there is a "
                                 "compositional field called peridotite."));

          if (this->convert_output_to_years() == true)
            melting_time_scale *= year_in_seconds;

          AssertThrow(melting_time_scale >= this->get_parameters().reaction_time_step,
                      ExcMessage("The reaction time step " + Utilities::to_string(this->get_parameters().reaction_time_step)
                                 + " in the operator splitting scheme is too large to compute melting rates! "
                                 "You have to choose it in such a way that it is smaller than the 'Melting time scale for "
                                 "operator splitting' chosen in the material model, which is currently "
                                 + Utilities::to_string(melting_time_scale) + "."));
          AssertThrow(melting_time_scale > 0,
                      ExcMessage("The Melting time scale for operator splitting must be larger than 0!"));

          if (this->include_melt_transport())
            {
              AssertThrow(this->introspection().compositional_name_exists("porosity"),
                          ExcMessage("Material model Melt visco plastic with melt transport only "
                                     "works if there is a compositional field called porosity."));
            }

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

    }

    template <int dim>
    void
    MeltViscoPlastic<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      rheology->create_plastic_outputs(out);

      if (rheology->use_elasticity)
        rheology->elastic_rheology.create_elastic_outputs(out);

      // if (out.template get_additional_output<PlasticAdditionalOutputs<dim> >() == nullptr)
      //   {
      //     const unsigned int n_points = out.n_evaluation_points();
      //     out.additional_outputs.push_back(
      //       std_cxx14::make_unique<MaterialModel::PlasticAdditionalOutputs<dim>> (n_points));
      //   }
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(MeltViscoPlastic,
                                   "EXPERIMENTAL melt visco plastic",
                                   "WARNING: This material model has not been thoroughly tested, "
                                   "and it is not advised that it be used for any serious models "
                                   "at this point.\n\n"
                                   "This material model implements a simple formulation of the "
                                   "material parameters required for the modelling of melt transport, "
                                   "including a source term for the porosity according to the melting "
                                   "model for dry peridotite of \\cite{KSL2003}. All other material "
                                   "properties are taken from the visco-plastic model.")
  }
}
