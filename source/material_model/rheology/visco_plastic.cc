/*
  Copyright (C) 2020 by the authors of the ASPECT code.

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

#include <aspect/material_model/rheology/visco_plastic.h>
#include <aspect/material_model/utilities.h>
#include <aspect/utilities.h>
#include <aspect/newton.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace
    {
      std::vector<std::string> make_plastic_additional_outputs_names()
      {
        std::vector<std::string> names;
        names.emplace_back("current_cohesions");
        names.emplace_back("current_friction_angles");
        names.emplace_back("plastic_yielding");
        return names;
      }
    }



    template <int dim>
    PlasticAdditionalOutputs<dim>::PlasticAdditionalOutputs (const unsigned int n_points)
      :
      NamedAdditionalMaterialOutputs<dim>(make_plastic_additional_outputs_names()),
      cohesions(n_points, numbers::signaling_nan<double>()),
      friction_angles(n_points, numbers::signaling_nan<double>()),
      yielding(n_points, numbers::signaling_nan<double>())
    {}



    template <int dim>
    std::vector<double>
    PlasticAdditionalOutputs<dim>::get_nth_output(const unsigned int idx) const
    {
      AssertIndexRange (idx, 3);
      switch (idx)
        {
          case 0:
            return cohesions;

          case 1:
            return friction_angles;

          case 2:
            return yielding;

          default:
            AssertThrow(false, ExcInternalError());
        }
      // We will never get here, so just return something
      return cohesions;
    }



    namespace Rheology
    {

      template <int dim>
      ViscoPlastic<dim>::ViscoPlastic ()
      {}



      template <int dim>
      std::pair<std::vector<double>, std::vector<bool> >
      ViscoPlastic<dim>::
      calculate_isostrain_viscosities (const MaterialModel::MaterialModelInputs<dim> &in,
                                       const unsigned int i,
                                       const std::vector<double> &volume_fractions,
                                       const std::vector<double> &phase_function_values,
                                       const std::vector<unsigned int> &n_phases_per_composition) const
      {

        // Initialize or fill variables used to calculate viscosities
        std::vector<bool> composition_yielding(volume_fractions.size(), false);
        std::vector<double> composition_viscosities(volume_fractions.size(), numbers::signaling_nan<double>());

        // Assemble stress tensor if elastic behavior is enabled
        SymmetricTensor<2,dim> stress_old = numbers::signaling_nan<SymmetricTensor<2,dim>>();
        if (use_elasticity == true)
          {
            for (unsigned int j=0; j < SymmetricTensor<2,dim>::n_independent_components; ++j)
              stress_old[SymmetricTensor<2,dim>::unrolled_to_component_indices(j)] = in.composition[i][j];
          }

        // The first time this function is called (first iteration of first time step)
        // a specified "reference" strain rate is used as the returned value would
        // otherwise be zero.
        const bool use_reference_strainrate = (this->get_timestep_number() == 0) &&
                                              (in.strain_rate[i].norm() <= std::numeric_limits<double>::min());

        double edot_ii;
        if (use_reference_strainrate)
          edot_ii = ref_strain_rate;
        else
          // Calculate the square root of the second moment invariant for the deviatoric strain rate tensor.
          edot_ii = std::max(std::sqrt(std::fabs(second_invariant(deviator(in.strain_rate[i])))),
                             min_strain_rate);

        // Calculate viscosities for each of the individual compositional phases
        for (unsigned int j=0; j < volume_fractions.size(); ++j)
          {
            // Step 1: viscous behavior

            // Choice of activation volume depends on whether there is an adiabatic temperature
            // gradient used when calculating the viscosity. This allows the same activation volume
            // to be used in incompressible and compressible models.
            const double temperature_for_viscosity = in.temperature[i] + adiabatic_temperature_gradient_for_viscosity*in.pressure[i];
            AssertThrow(temperature_for_viscosity != 0, ExcMessage(
                          "The temperature used in the calculation of the visco-plastic rheology is zero. "
                          "This is not allowed, because this value is used to divide through. It is probably "
                          "being caused by the temperature being zero somewhere in the model. The relevant "
                          "values for debugging are: temperature (" + Utilities::to_string(in.temperature[i]) +
                          "), adiabatic_temperature_gradient_for_viscosity ("
                          + Utilities::to_string(adiabatic_temperature_gradient_for_viscosity) + ") and pressure ("
                          + Utilities::to_string(in.pressure[i]) + ")."));

            // Step 1a: compute viscosity from diffusion creep law
            const double viscosity_diffusion = diffusion_creep.compute_viscosity(in.pressure[i], temperature_for_viscosity, j,
                                                                                 phase_function_values,
                                                                                 n_phases_per_composition);

            // Step 1b: compute viscosity from dislocation creep law
            const double viscosity_dislocation = dislocation_creep.compute_viscosity(edot_ii, in.pressure[i], temperature_for_viscosity, j,
                                                                                     phase_function_values,
                                                                                     n_phases_per_composition);

            // Step 1c: select what form of viscosity to use (diffusion, dislocation, fk, or composite)
            double viscosity_pre_yield = 0.0;
            switch (viscous_flow_law)
              {
                case diffusion:
                {
                  viscosity_pre_yield = viscosity_diffusion;
                  break;
                }
                case dislocation:
                {
                  viscosity_pre_yield = viscosity_dislocation;
                  break;
                }
                case frank_kamenetskii:
                {
                  viscosity_pre_yield = frank_kamenetskii_rheology->compute_viscosity(in.temperature[i], j);
                  break;
                }
                case composite:
                {
                  viscosity_pre_yield = (viscosity_diffusion * viscosity_dislocation)/
                                        (viscosity_diffusion + viscosity_dislocation);
                  break;
                }
                default:
                {
                  AssertThrow(false, ExcNotImplemented());
                  break;
                }
              }

            // Step 1d: compute viscosity from Peierls creep law and harmonically average with current viscosities
            if (use_peierls_creep)
              {
                const double viscosity_peierls = peierls_creep->compute_viscosity(edot_ii, in.pressure[i], temperature_for_viscosity, j,
                                                                                  phase_function_values,
                                                                                  n_phases_per_composition);
                viscosity_pre_yield = (viscosity_pre_yield * viscosity_peierls) / (viscosity_pre_yield + viscosity_peierls);
              }

            // Step 1e: multiply the viscosity by a constant (default value is 1)
            viscosity_pre_yield = constant_viscosity_prefactors.compute_viscosity(viscosity_pre_yield, j);

            // Step 2: calculate the viscous stress magnitude
            // and strain rate. If requested compute visco-elastic contributions.
            double current_edot_ii = edot_ii;

            if (use_elasticity)
              {
                const std::vector<double> &elastic_shear_moduli = elastic_rheology.get_elastic_shear_moduli();

                if (use_reference_strainrate == true)
                  current_edot_ii = ref_strain_rate;
                else
                  {
                    const double viscoelastic_strain_rate_invariant = elastic_rheology.calculate_viscoelastic_strain_rate(in.strain_rate[i],
                                                                      stress_old,
                                                                      elastic_shear_moduli[j]);

                    current_edot_ii = std::max(viscoelastic_strain_rate_invariant,
                                               min_strain_rate);
                  }

                // Step 2a: calculate viscoelastic (effective) viscosity
                viscosity_pre_yield = elastic_rheology.calculate_viscoelastic_viscosity(viscosity_pre_yield,
                                                                                        elastic_shear_moduli[j]);
              }

            // Step 2b: calculate current (viscous or viscous + elastic) stress magnitude
            double current_stress = 2. * viscosity_pre_yield * current_edot_ii;

            // Step 3: strain weakening

            // Step 3a: calculate strain weakening factors for the cohesion, friction, and pre-yield viscosity
            // If no brittle and/or viscous strain weakening is applied, the factors are 1.
            const std::array<double, 3> weakening_factors = strain_rheology.compute_strain_weakening_factors(j, in.composition[i]);

            // Step 3b: calculate weakened friction, cohesion, and pre-yield viscosity and adjust the current_stress accordingly
            const double current_cohesion = drucker_prager_parameters.cohesions[j] * weakening_factors[0];
            const double current_friction = drucker_prager_parameters.angles_internal_friction[j] * weakening_factors[1];
            viscosity_pre_yield *= weakening_factors[2];
            current_stress *= weakening_factors[2];

            // Step 4: plastic yielding

            // Determine if the pressure used in Drucker Prager plasticity will be capped at 0 (default).
            // This may be necessary in models without gravity and the dynamic stresses are much higher
            // than the lithostatic pressure.

            double pressure_for_plasticity = in.pressure[i];
            if (allow_negative_pressures_in_plasticity == false)
              pressure_for_plasticity = std::max(in.pressure[i],0.0);

            // Step 4a: calculate Drucker-Prager yield stress
            const double yield_stress = drucker_prager_plasticity.compute_yield_stress(current_cohesion,
                                                                                       current_friction,
                                                                                       pressure_for_plasticity,
                                                                                       drucker_prager_parameters.max_yield_stress);

            // Step 4b: select if yield viscosity is based on Drucker Prager or stress limiter rheology
            double viscosity_yield = viscosity_pre_yield;
            switch (yield_mechanism)
              {
                case stress_limiter:
                {
                  //Step 4b-1: always rescale the viscosity back to the yield surface
                  const double viscosity_limiter = yield_stress / (2.0 * ref_strain_rate)
                                                   * std::pow((edot_ii/ref_strain_rate),
                                                              1./exponents_stress_limiter[j] - 1.0);
                  viscosity_yield = 1. / ( 1./viscosity_limiter + 1./viscosity_pre_yield);
                  break;
                }
                case drucker_prager:
                {
                  // Step 4b-2: if the current stress is greater than the yield stress,
                  // rescale the viscosity back to yield surface
                  if (current_stress >= yield_stress)
                    {
                      viscosity_yield = drucker_prager_plasticity.compute_viscosity(current_cohesion,
                                                                                    current_friction,
                                                                                    pressure_for_plasticity,
                                                                                    current_edot_ii,
                                                                                    drucker_prager_parameters.max_yield_stress,
                                                                                    viscosity_pre_yield);
                      composition_yielding[j] = true;
                    }
                  break;
                }
                default:
                {
                  AssertThrow(false, ExcNotImplemented());
                  break;
                }
              }

            // Step 5: limit the viscosity with specified minimum and maximum bounds
            composition_viscosities[j] = std::min(std::max(viscosity_yield, min_visc), max_visc);
          }
        return std::make_pair (composition_viscosities, composition_yielding);
      }



      template <int dim>
      void
      ViscoPlastic<dim>::
      compute_viscosity_derivatives(const unsigned int i,
                                    const std::vector<double> &volume_fractions,
                                    const std::vector<double> &composition_viscosities,
                                    const MaterialModel::MaterialModelInputs<dim> &in,
                                    MaterialModel::MaterialModelOutputs<dim> &out,
                                    const std::vector<double> &phase_function_values,
                                    const std::vector<unsigned int> &n_phases_per_composition) const
      {
        MaterialModel::MaterialModelDerivatives<dim> *derivatives =
          out.template get_additional_output<MaterialModel::MaterialModelDerivatives<dim> >();

        if (derivatives != nullptr)
          {
            // compute derivatives if necessary
            std::vector<SymmetricTensor<2,dim> > composition_viscosities_derivatives(volume_fractions.size());
            std::vector<double> composition_dviscosities_dpressure(volume_fractions.size());

            const double finite_difference_accuracy = 1e-7;

            // A new material model inputs variable that uses the strain rate and pressure difference.
            MaterialModel::MaterialModelInputs<dim> in_derivatives = in;

            // For each independent component, compute the derivative.
            for (unsigned int component = 0; component < SymmetricTensor<2,dim>::n_independent_components; ++component)
              {
                const TableIndices<2> strain_rate_indices = SymmetricTensor<2,dim>::unrolled_to_component_indices (component);

                // components that are not on the diagonal are multiplied by 0.5, because the symmetric tensor
                // is modified by 0.5 in both symmetric directions (xy/yx) simultaneously and we compute the combined
                // derivative
                const SymmetricTensor<2,dim> strain_rate_difference = in.strain_rate[i]
                                                                      + std::max(std::fabs(in.strain_rate[i][strain_rate_indices]), min_strain_rate)
                                                                      * (component > dim-1 ? 0.5 : 1 )
                                                                      * finite_difference_accuracy
                                                                      * Utilities::nth_basis_for_symmetric_tensors<dim>(component);

                in_derivatives.strain_rate[i] = strain_rate_difference;

                std::vector<double> eta_component =
                  calculate_isostrain_viscosities(in_derivatives, i, volume_fractions,
                                                  phase_function_values, n_phases_per_composition).first;

                // For each composition of the independent component, compute the derivative.
                for (unsigned int composition_index = 0; composition_index < eta_component.size(); ++composition_index)
                  {
                    // compute the difference between the viscosity with and without the strain-rate difference.
                    double viscosity_derivative = eta_component[composition_index] - composition_viscosities[composition_index];
                    if (viscosity_derivative != 0)
                      {
                        // when the difference is non-zero, divide by the difference.
                        viscosity_derivative /= std::max(std::fabs(strain_rate_difference[strain_rate_indices]), min_strain_rate)
                                                * finite_difference_accuracy;
                      }
                    composition_viscosities_derivatives[composition_index][strain_rate_indices] = viscosity_derivative;
                  }
              }

            /**
             * Now compute the derivative of the viscosity to the pressure
             */
            const double pressure_difference = in.pressure[i] + (std::fabs(in.pressure[i]) * finite_difference_accuracy);

            in_derivatives.pressure[i] = pressure_difference;

            // Modify the in_derivatives object again to take the original strain rate.
            in_derivatives.strain_rate[i] = in.strain_rate[i];

            const std::vector<double> viscosity_difference =
              calculate_isostrain_viscosities(in_derivatives, i, volume_fractions,
                                              phase_function_values, n_phases_per_composition).first;

            for (unsigned int composition_index = 0; composition_index < viscosity_difference.size(); ++composition_index)
              {
                double viscosity_derivative = viscosity_difference[composition_index] - composition_viscosities[composition_index];
                if (viscosity_difference[composition_index] != 0)
                  {
                    if (in.pressure[i] != 0)
                      {
                        viscosity_derivative /= std::fabs(in.pressure[i]) * finite_difference_accuracy;
                      }
                    else
                      {
                        viscosity_derivative = 0;
                      }
                  }
                composition_dviscosities_dpressure[composition_index] = viscosity_derivative;
              }

            double viscosity_averaging_p = 0; // Geometric
            if (viscosity_averaging == MaterialUtilities::harmonic)
              viscosity_averaging_p = -1;
            if (viscosity_averaging == MaterialUtilities::arithmetic)
              viscosity_averaging_p = 1;
            if (viscosity_averaging == MaterialUtilities::maximum_composition)
              viscosity_averaging_p = 1000;

            derivatives->viscosity_derivative_wrt_strain_rate[i] =
              Utilities::derivative_of_weighted_p_norm_average(out.viscosities[i],
                                                               volume_fractions,
                                                               composition_viscosities,
                                                               composition_viscosities_derivatives,
                                                               viscosity_averaging_p);
            derivatives->viscosity_derivative_wrt_pressure[i] =
              Utilities::derivative_of_weighted_p_norm_average(out.viscosities[i],
                                                               volume_fractions,
                                                               composition_viscosities,
                                                               composition_dviscosities_dpressure,
                                                               viscosity_averaging_p);
          }
      }



      template <int dim>
      ComponentMask
      ViscoPlastic<dim>::
      get_volumetric_composition_mask() const
      {
        // Store which components to exclude during the volume fraction computation.
        ComponentMask composition_mask = strain_rheology.get_strain_composition_mask();

        if (use_elasticity)
          {
            for (unsigned int i = 0; i < SymmetricTensor<2,dim>::n_independent_components ; ++i)
              composition_mask.set(i,false);
          }

        return composition_mask;
      }



      template <int dim>
      void
      ViscoPlastic<dim>::declare_parameters (ParameterHandler &prm)
      {
        Rheology::StrainDependent<dim>::declare_parameters (prm);

        Rheology::Elasticity<dim>::declare_parameters (prm);

        // Reference and minimum/maximum values
        prm.declare_entry ("Minimum strain rate", "1.0e-20", Patterns::Double (0.),
                           "Stabilizes strain dependent viscosity. Units: \\si{\\per\\second}.");
        prm.declare_entry ("Reference strain rate","1.0e-15",Patterns::Double (0.),
                           "Reference strain rate for first time step. Units: \\si{\\per\\second}.");
        prm.declare_entry ("Minimum viscosity", "1e17", Patterns::Double (0.),
                           "Lower cutoff for effective viscosity. Units: \\si{\\pascal\\second}.");
        prm.declare_entry ("Maximum viscosity", "1e28", Patterns::Double (0.),
                           "Upper cutoff for effective viscosity. Units: \\si{\\pascal\\second}.");
        prm.declare_entry ("Reference viscosity", "1e22", Patterns::Double (0.),
                           "Reference viscosity for nondimensionalization. "
                           "To understand how pressure scaling works, take a look at "
                           "\\cite{KHB12}. In particular, the value of this parameter "
                           "would not affect the solution computed by \\aspect{} if "
                           "we could do arithmetic exactly; however, computers do "
                           "arithmetic in finite precision, and consequently we need to "
                           "scale quantities in ways so that their magnitudes are "
                           "roughly the same. As explained in \\cite{KHB12}, we scale "
                           "the pressure during some computations (never visible by "
                           "users) by a factor that involves a reference viscosity. This "
                           "parameter describes this reference viscosity."
                           "\n\n"
                           "For problems with a constant viscosity, you will generally want "
                           "to choose the reference viscosity equal to the actual viscosity. "
                           "For problems with a variable viscosity, the reference viscosity "
                           "should be a value that adequately represents the order of "
                           "magnitude of the viscosities that appear, such as an average "
                           "value or the value one would use to compute a Rayleigh number."
                           "\n\n"
                           "Units: \\si{\\pascal\\second}.");

        // Rheological parameters
        prm.declare_entry ("Viscosity averaging scheme", "harmonic",
                           Patterns::Selection("arithmetic|harmonic|geometric|maximum composition"),
                           "When more than one compositional field is present at a point "
                           "with different viscosities, we need to come up with an average "
                           "viscosity at that point.  Select a weighted harmonic, arithmetic, "
                           "geometric, or maximum composition.");
        prm.declare_entry ("Viscous flow law", "composite",
                           Patterns::Selection("diffusion|dislocation|frank kamenetskii|composite"),
                           "Select what type of viscosity law to use between diffusion, "
                           "dislocation, frank kamenetskii, and composite options. Soon there will be an option "
                           "to select a specific flow law for each assigned composition ");
        prm.declare_entry ("Yield mechanism", "drucker",
                           Patterns::Selection("drucker|limiter"),
                           "Select what type of yield mechanism to use between Drucker Prager "
                           "and stress limiter options.");
        prm.declare_entry ("Allow negative pressures in plasticity", "false",
                           Patterns::Bool (),
                           "Whether to allow negative pressures to be used in the computation "
                           "of plastic yield stresses and viscosities. Setting this parameter "
                           "to true may be advantageous in models without gravity where the "
                           "dynamic stresses are much higher than the lithostatic pressure. "
                           "If false, the minimum pressure in the plasticity formulation will "
                           "be set to zero.");

        // Diffusion creep parameters
        Rheology::DiffusionCreep<dim>::declare_parameters(prm);

        // Dislocation creep parameters
        Rheology::DislocationCreep<dim>::declare_parameters(prm);

        // Frank-Kamenetskii viscosity parameters
        Rheology::FrankKamenetskii<dim>::declare_parameters(prm);

        // Peierls creep parameters
        Rheology::PeierlsCreep<dim>::declare_parameters(prm);

        prm.declare_entry ("Include Peierls creep", "false",
                           Patterns::Bool (),
                           "Whether to include Peierls creep in the rheological formulation.");

        // Constant viscosity prefactor parameters
        Rheology::ConstantViscosityPrefactors<dim>::declare_parameters(prm);

        // Drucker Prager plasticity parameters
        Rheology::DruckerPrager<dim>::declare_parameters(prm);

        // Stress limiter parameters
        prm.declare_entry ("Stress limiter exponents", "1.0",
                           Patterns::List(Patterns::Double (0.)),
                           "List of stress limiter exponents, $n_{\\text{lim}}$, "
                           "for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "Units: none.");

        // Temperature gradient in viscosity laws to include an adiabat (note units of K/Pa)
        prm.declare_entry ("Adiabat temperature gradient for viscosity", "0.0", Patterns::Double (0.),
                           "Add an adiabatic temperature gradient to the temperature used in the flow law "
                           "so that the activation volume is consistent with what one would use in a "
                           "earth-like (compressible) model. Default is set so this is off. "
                           "Note that this is a linear approximation of the real adiabatic gradient, which "
                           "is okay for the upper mantle, but is not really accurate for the lower mantle. "
                           "Using a pressure gradient of 32436 Pa/m, then a value of "
                           "0.3 K/km = 0.0003 K/m = 9.24e-09 K/Pa gives an earth-like adiabat."
                           "Units: \\si{\\kelvin\\per\\pascal}.");

        prm.declare_entry ("Include viscoelasticity", "false",
                           Patterns::Bool (),
                           "Whether to include elastic effects in the rheological formulation.");
      }



      template <int dim>
      void
      ViscoPlastic<dim>::parse_parameters (ParameterHandler &prm,
                                           const std::shared_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition)
      {

        // increment by one for background:
        const unsigned int n_fields = this->n_compositional_fields() + 1;

        strain_rheology.initialize_simulator (this->get_simulator());
        strain_rheology.parse_parameters(prm);

        use_elasticity = prm.get_bool ("Include viscoelasticity");

        if (use_elasticity)
          {
            elastic_rheology.initialize_simulator (this->get_simulator());
            elastic_rheology.parse_parameters(prm);
          }

        // Reference and minimum/maximum values
        min_strain_rate = prm.get_double("Minimum strain rate");
        ref_strain_rate = prm.get_double("Reference strain rate");
        min_visc = prm.get_double ("Minimum viscosity");
        max_visc = prm.get_double ("Maximum viscosity");
        ref_visc = prm.get_double ("Reference viscosity");

        viscosity_averaging = MaterialUtilities::parse_compositional_averaging_operation ("Viscosity averaging scheme",
                              prm);

        // Rheological parameters
        if (prm.get ("Viscous flow law") == "composite")
          viscous_flow_law = composite;
        else if (prm.get ("Viscous flow law") == "diffusion")
          viscous_flow_law = diffusion;
        else if (prm.get ("Viscous flow law") == "dislocation")
          viscous_flow_law = dislocation;
        else if (prm.get ("Viscous flow law") == "frank kamenetskii")
          viscous_flow_law = frank_kamenetskii;
        else
          AssertThrow(false, ExcMessage("Not a valid viscous flow law"));

        if (prm.get ("Yield mechanism") == "drucker")
          yield_mechanism = drucker_prager;
        else if (prm.get ("Yield mechanism") == "limiter")
          yield_mechanism = stress_limiter;
        else
          AssertThrow(false, ExcMessage("Not a valid yield mechanism."));

        AssertThrow(use_elasticity == false || yield_mechanism == drucker_prager,
                    ExcMessage("Elastic behavior is only tested with the "
                               "'drucker prager' plasticity option."));

        allow_negative_pressures_in_plasticity = prm.get_bool ("Allow negative pressures in plasticity");

        // Diffusion creep parameters
        diffusion_creep.initialize_simulator (this->get_simulator());
        diffusion_creep.parse_parameters(prm, expected_n_phases_per_composition);

        // Dislocation creep parameters
        dislocation_creep.initialize_simulator (this->get_simulator());
        dislocation_creep.parse_parameters(prm, expected_n_phases_per_composition);

        // Frank Kamenetskii viscosity parameters
        if (viscous_flow_law == frank_kamenetskii)
          {
            frank_kamenetskii_rheology = std_cxx14::make_unique<Rheology::FrankKamenetskii<dim>>();
            frank_kamenetskii_rheology->initialize_simulator (this->get_simulator());
            frank_kamenetskii_rheology->parse_parameters(prm);
          }

        // Peierls creep parameters
        use_peierls_creep = prm.get_bool ("Include Peierls creep");
        if (use_peierls_creep)
          {
            peierls_creep = std_cxx14::make_unique<Rheology::PeierlsCreep<dim>>();
            peierls_creep->initialize_simulator (this->get_simulator());
            peierls_creep->parse_parameters(prm);
          }

        // Constant viscosity prefactor parameters
        constant_viscosity_prefactors.initialize_simulator (this->get_simulator());
        constant_viscosity_prefactors.parse_parameters(prm);

        // Plasticity parameters
        drucker_prager_parameters = drucker_prager_plasticity.parse_parameters(this->n_compositional_fields()+1,
                                                                               prm);

        // Stress limiter parameter
        exponents_stress_limiter  = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Stress limiter exponents"))),
                                                                            n_fields,
                                                                            "Stress limiter exponents");

        // Include an adiabat temperature gradient in flow laws
        adiabatic_temperature_gradient_for_viscosity = prm.get_double("Adiabat temperature gradient for viscosity");
        if (this->get_heating_model_manager().adiabatic_heating_enabled())
          AssertThrow (adiabatic_temperature_gradient_for_viscosity == 0.0,
                       ExcMessage("If adiabatic heating is enabled you should not add another adiabatic gradient"
                                  "to the temperature for computing the viscosity, because the ambient"
                                  "temperature profile already includes the adiabatic gradient."));

      }



      template <int dim>
      void
      ViscoPlastic<dim>::create_plastic_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        if (out.template get_additional_output<PlasticAdditionalOutputs<dim> >() == nullptr)
          {
            const unsigned int n_points = out.n_evaluation_points();
            out.additional_outputs.push_back(
              std_cxx14::make_unique<PlasticAdditionalOutputs<dim>> (n_points));
          }
      }



      template <int dim>
      void
      ViscoPlastic<dim>::
      fill_plastic_outputs(const unsigned int i,
                           const std::vector<double> &volume_fractions,
                           const bool plastic_yielding,
                           const MaterialModel::MaterialModelInputs<dim> &in,
                           MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        PlasticAdditionalOutputs<dim> *plastic_out = out.template get_additional_output<PlasticAdditionalOutputs<dim> >();

        if (plastic_out != nullptr)
          {
            plastic_out->cohesions[i] = 0;
            plastic_out->friction_angles[i] = 0;
            plastic_out->yielding[i] = plastic_yielding ? 1 : 0;

            // set to weakened values, or unweakened values when strain weakening is not used
            for (unsigned int j=0; j < volume_fractions.size(); ++j)
              {
                // Calculate the strain weakening factors and weakened values
                const std::array<double, 3> weakening_factors = strain_rheology.compute_strain_weakening_factors(j, in.composition[i]);
                plastic_out->cohesions[i]   += volume_fractions[j] * (drucker_prager_parameters.cohesions[j] * weakening_factors[0]);
                // Also convert radians to degrees
                plastic_out->friction_angles[i] += 180.0/numbers::PI * volume_fractions[j] * (drucker_prager_parameters.angles_internal_friction[j] * weakening_factors[1]);
              }
          }
      }

    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  namespace Rheology \
  { \
    template class ViscoPlastic<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
