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


#include <aspect/material_model/reaction_model/grain_size_evolution.h>
#include <aspect/utilities.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/heating_model/shear_heating.h>
#include <aspect/simulator_signals.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/sundials/arkode.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      // namespace for helper functions
      namespace
      {
        // Compute the nth moment of the log-normal distribution as used in
        // Bercovici and Ricard (2012; F6) using mean = 0 and variance = 0.8.
        double nth_moment_of_lognormal_distribution (const unsigned int n)
        {
          const double sigma = 0.8;
          return std::exp(n * n * sigma * sigma / 2.);
        }



        // Computes the product of two phase fractions (assuming exactly two phases)
        // as used in Bercovici and Ricard (2012).
        double phase_distribution_function (const double volume_fraction_phase_one)
        {
          const double volume_fraction_phase_two = 1. - volume_fraction_phase_one;
          return (volume_fraction_phase_one * volume_fraction_phase_two);
        }



        // Compute the factor to convert from an interface roughness equation
        // to a mean grain size. See Appendix H.1, Eqs. 8 and F.28 in
        // Bercovici and Ricard (2012) for more details.
        double
        roughness_to_grain_size_factor (const double volume_fraction_phase_one)
        {
          const double b1 = 1./20 ;
          const double c1 = 3.0 * b1 * nth_moment_of_lognormal_distribution(4) / (8.0 * nth_moment_of_lognormal_distribution (2));

          const double volume_fraction_phase_two = 1. - volume_fraction_phase_one;

          const double h1 = c1 * (1 - volume_fraction_phase_one);
          const double h2 = c1 * (1 - volume_fraction_phase_two);

          const double one_over_sqrt_h = volume_fraction_phase_one / std::sqrt(h1) + volume_fraction_phase_two / std::sqrt(h2);

          return (1./one_over_sqrt_h);
        }
      }



      template <int dim>
      void
      GrainSizeEvolution<dim>::initialize_phase_function(std::shared_ptr<MaterialUtilities::PhaseFunction<dim>> &phase_function_)
      {
        CitationInfo::add("grainsize");

        phase_function = phase_function_;
        n_phase_transitions = phase_function->n_phases_for_each_composition()[0] - 1;
      }



      template <int dim>
      double
      GrainSizeEvolution<dim>::
      compute_partitioning_fraction (const double temperature) const
      {
        const double power_term_base = maximum_grain_size_reduction_work_fraction/minimum_grain_size_reduction_work_fraction;

        const double power_term_numerator    =  temperature_minimum_partitioning_power -
                                                std::pow (temperature, grain_size_reduction_work_fraction_exponent);

        const double power_term_denominator  =  temperature_minimum_partitioning_power -
                                                temperature_maximum_partitioning_power;

        // We have to ensure the power term exponent is between 0 and 1, otherwise the partitioning fraction
        // will be outside the set bounds for the work fraction.
        const double power_term_exponent = std::max(std::min(power_term_numerator / power_term_denominator, 1.0), 0.0);

        const double power_term = std::pow(power_term_base,
                                           power_term_exponent);

        return minimum_grain_size_reduction_work_fraction * power_term;
      }



      template <int dim>
      void
      GrainSizeEvolution<dim>::
      calculate_reaction_terms (const typename Interface<dim>::MaterialModelInputs  &in,
                                const std::vector<double>                           &pressures,
                                const std::vector<unsigned int>                     &phase_indices,
                                const std::function<double(const double, const double,
                                                           const double,const SymmetricTensor<2,dim> &,const unsigned int,const double,const double)>          &dislocation_viscosity,
                                const std::function<double(
                                  const double,const double,const double,const double,const double,const unsigned int)> &diffusion_viscosity,
                                const double                                         min_eta,
                                const double                                         max_eta,
                                typename Interface<dim>::MaterialModelOutputs       &out) const
      {
        // We want to iterate over the grain size evolution here, as we solve in fact an ordinary differential equation
        // and it is not correct to use the starting grain size (which also introduces instabilities).
        // We assume that the strain rate is constant across all substeps for the ODE (within one advection step).
        // Even though this assumption may not always be the most accurate (compared to assuming a constant stress),
        // it leads to a more stable behavior because it implies that a reduction in grain size leads to less work
        // being done by dislocation creep, so that the grain size reduction rate is lower in the next substep.
        std::vector<std::vector<double>> &reaction_terms = out.reaction_terms;

        const unsigned int n_evaluation_points = in.n_evaluation_points();
        const unsigned int grain_size_index = this->introspection().get_indices_for_fields_of_type(CompositionalFieldDescription::grain_size)[0];

        // Set up a vector that tells us which phase transition has been crossed for each point we are evaluating.
        std::vector<int> crossed_transitions (n_evaluation_points, -1);

        // Copy grain sizes into a VectorType for the ODE solver.
        // While there also limit grain sizes to be larger than the minimum grain size.
        using VectorType = Vector<double>;
        VectorType grain_sizes(n_evaluation_points);
        for (unsigned int i=0; i<n_evaluation_points; ++i)
          grain_sizes[i] = std::max(minimum_grain_size, in.composition[i][grain_size_index]);

        const double timestep = this->simulator_is_past_initialization()
                                ?
                                this->get_timestep()
                                :
                                0.0;

        // If the grain sizes are not valid, we are in initialization, or the time
        // step size is zero, we do not solve the ODE.
        if (std::all_of(grain_sizes.begin(), grain_sizes.end(), [](double gs)
        {
          return gs != gs || gs < std::numeric_limits<double>::min();
          })
        || timestep == 0.0)
        return;

        SUNDIALS::ARKode<VectorType>::AdditionalData data;

        data.initial_time = 0.0;
        data.final_time = this->get_timestep();

        // TODO: make this an input parameter.
        data.initial_step_size = 0.001 * this->get_timestep();
        data.output_period = this->get_timestep();
        data.minimum_step_size = 1.e-6*this->get_timestep();
        data.maximum_order = 3;
        data.maximum_non_linear_iterations = 30;

        // Because both tolerances are added, we set the absolute
        // tolerance to 0.
        data.relative_tolerance = 1e-3;
        data.absolute_tolerance = 0;

        SUNDIALS::ARKode<VectorType> ode(data);
        ode.explicit_function = [&] (const double     /*time*/,
                                     const VectorType &y,
                                     VectorType       &grain_size_rates_of_change)
        {
          for (unsigned int i=0; i<n_evaluation_points; ++i)
            {
              const double grain_size = std::max(minimum_grain_size, y[i]);

              // Precompute the partitioning_fraction since it is constant during the evolution.
              // This is only used for the pinned_grain_damage formulation.
              const double partitioning_fraction = (grain_size_evolution_formulation == Formulation::pinned_grain_damage)
                                                   ?
                                                   compute_partitioning_fraction(in.temperature[i])
                                                   :
                                                   0.0;

              // We keep the dislocation viscosity of the last iteration as guess
              // for the next one.
              double current_dislocation_viscosity = 0.0;

              const double adiabatic_temperature = this->get_adiabatic_conditions().is_initialized()
                                                   ?
                                                   this->get_adiabatic_conditions().temperature(in.position[i])
                                                   :
                                                   in.temperature[i];



              // grain size growth due to Ostwald ripening
              const double m = grain_growth_exponent[phase_indices[i]];

              double grain_size_growth_rate = grain_growth_rate_constant[phase_indices[i]] / (m * std::pow(grain_size,m-1))
                                              * std::exp(- (grain_growth_activation_energy[phase_indices[i]] + pressures[i] * grain_growth_activation_volume[phase_indices[i]])
                                                         / (constants::gas_constant * in.temperature[i]));

              // in the two-phase damage model grain growth depends on the proportion of the two phases
              if (grain_size_evolution_formulation == Formulation::pinned_grain_damage)
                grain_size_growth_rate *= geometric_constant[phase_indices[i]] * phase_distribution /
                                          std::pow(roughness_to_grain_size, m);

              // grain size reduction in dislocation creep regime
              const SymmetricTensor<2,dim> shear_strain_rate = in.strain_rate[i] - 1./dim * trace(in.strain_rate[i]) * unit_symmetric_tensor<dim>();
              const double second_strain_rate_invariant = std::sqrt(std::max(-second_invariant(shear_strain_rate), 0.));

              const double current_diffusion_viscosity   = diffusion_viscosity(in.temperature[i], adiabatic_temperature, pressures[i], grain_size, second_strain_rate_invariant, phase_indices[i]);
              current_dislocation_viscosity = dislocation_viscosity(in.temperature[i], adiabatic_temperature, pressures[i], in.strain_rate[i], phase_indices[i], current_diffusion_viscosity, current_dislocation_viscosity);

              double current_viscosity;
              if (std::abs(second_strain_rate_invariant) > 1e-30)
                current_viscosity = current_dislocation_viscosity * current_diffusion_viscosity / (current_dislocation_viscosity + current_diffusion_viscosity);
              else
                current_viscosity = current_diffusion_viscosity;

              const double dislocation_strain_rate = second_strain_rate_invariant
                                                     * current_viscosity / current_dislocation_viscosity;

              double grain_size_reduction_rate = 0.0;

              if (grain_size_evolution_formulation == Formulation::paleowattmeter)
                {
                  // paleowattmeter: Austin and Evans (2007): Paleowattmeters: A scaling relation for dynamically recrystallized grain size. Geology 35, 343-346
                  const double stress = 2.0 * second_strain_rate_invariant * std::min(std::max(min_eta,current_viscosity),max_eta);
                  grain_size_reduction_rate = 2.0 * stress * boundary_area_change_work_fraction[phase_indices[i]] * dislocation_strain_rate * grain_size * grain_size
                                              / (geometric_constant[phase_indices[i]] * grain_boundary_energy[phase_indices[i]]);
                }
              else if (grain_size_evolution_formulation == Formulation::pinned_grain_damage)
                {
                  // pinned_grain_damage: Mulyukova and Bercovici (2018) Collapse of passive margins by lithospheric damage and plunging grain size. Earth and Planetary Science Letters, 484, 341-352.
                  const double stress = 2.0 * second_strain_rate_invariant * std::min(std::max(min_eta,current_viscosity),max_eta);
                  grain_size_reduction_rate = 2.0 * stress * partitioning_fraction * second_strain_rate_invariant * grain_size * grain_size
                                              * roughness_to_grain_size
                                              / (geometric_constant[phase_indices[i]] * grain_boundary_energy[phase_indices[i]] * phase_distribution);
                }
              else if (grain_size_evolution_formulation == Formulation::paleopiezometer)
                {
                  // paleopiezometer: Hall and Parmentier (2003): Influence of grain size evolution on convective instability. Geochem. Geophys. Geosyst., 4(3).
                  grain_size_reduction_rate = reciprocal_required_strain[phase_indices[i]] * dislocation_strain_rate * grain_size;
                }
              else
                AssertThrow(false, ExcNotImplemented());

              grain_size_rates_of_change[i] = grain_size_growth_rate - grain_size_reduction_rate;
            }
        };


        const unsigned int iteration_count = ode.solve_ode(grain_sizes);
        this->get_signals().post_ARKode_solve(*this, iteration_count);

        for (unsigned int i=0; i<n_evaluation_points; ++i)
          {
            Assert(grain_sizes[i] > 0,
                   ExcMessage("The grain size became smaller than zero. This is not valid, "
                              "and likely an effect of a too large sub-timestep, or unrealistic "
                              "input parameters."));


            // Reduce grain size to recrystallized_grain_size if the grain crossed a phase transition.
            // To do so, check if a grain has moved further than its distance to a phase transition in
            // vertical direction and the movement is away from the phase transition.
            // 'crossed_transition' will be -1 if we crossed no transition, or the index of a
            // phase transition, if we crossed one.
            int crossed_transition = -1;

            const double gravity_norm = this->get_gravity_model().gravity_vector(in.position[i]).norm();
            Tensor<1,dim> vertical_direction = this->get_gravity_model().gravity_vector(in.position[i]);
            if (gravity_norm > 0.0)
              vertical_direction /= gravity_norm;

            for (unsigned int phase=0; phase<n_phase_transitions; ++phase)
              {
                // Both distances are positive when they are downward from the transition (since gravity points down)
                const double distance_from_transition = this->get_geometry_model().depth(in.position[i]) - phase_function->get_transition_depth(phase);
                const double distance_moved = in.velocity[i] * vertical_direction * timestep;

                // To make sure we actually reset the grain size of all the material passing through
                // the transition, we take 110% of the distance a grain has moved for the check.
                if (std::abs(distance_moved) * 1.1 > std::abs(distance_from_transition)
                    &&
                    distance_moved * distance_from_transition >= 0)
                  crossed_transition = phase;
              }

            // TODO: recrystallize at correct time while doing grain size evolution instead of afterwards
            double phase_grain_size_reduction = 0.0;
            if (this->get_timestep_number() > 0)
              {
                // check if material has crossed any phase transition, if yes, reset grain size
                if (crossed_transition != -1)
                  if (recrystallized_grain_size[crossed_transition] > 0.0)
                    phase_grain_size_reduction = grain_sizes[i] - recrystallized_grain_size[crossed_transition];
              }

            grain_sizes[i] = std::max(grain_sizes[i], minimum_grain_size);
            const double grain_size_change = grain_sizes[i] - in.composition[i][grain_size_index] - phase_grain_size_reduction;

            // this reaction model is only responsible for the grain size field
            reaction_terms[i][grain_size_index] = grain_size_change;
          }

        return;
      }




      template <int dim>
      void
      GrainSizeEvolution<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Grain growth activation energy", "3.5e5",
                           Patterns::List (Patterns::Double (0.)),
                           "The activation energy for grain growth $E_g$. "
                           "List must have one more entry than the Phase transition depths. "
                           "Units: \\si{\\joule\\per\\mole}.");
        prm.declare_entry ("Grain growth activation volume", "8e-6",
                           Patterns::List (Patterns::Double (0.)),
                           "The activation volume for grain growth $V_g$. "
                           "List must have one more entry than the Phase transition depths. "
                           "Units: \\si{\\meter\\cubed\\per\\mole}.");
        prm.declare_entry ("Grain growth exponent", "3.",
                           Patterns::List (Patterns::Double (0.)),
                           "The exponent of the grain growth law $p_g$. This is an experimentally determined "
                           "grain growth constant. "
                           "List must have one more entry than the Phase transition depths. "
                           "Units: none.");
        prm.declare_entry ("Grain growth rate constant", "1.5e-5",
                           Patterns::List (Patterns::Double (0.)),
                           "The prefactor for the Ostwald ripening grain growth law $G_0$. "
                           "This is dependent on water content, which is assumed to be "
                           "50 H/$10^6$ Si for the default value. "
                           "List must have one more entry than the Phase transition depths. "
                           "Units: \\si{\\meter}$^{p_g}$\\si{\\per\\second}.");
        prm.declare_entry ("Reciprocal required strain", "10.",
                           Patterns::List (Patterns::Double (0.)),
                           "This parameter ($\\lambda$) gives an estimate of the strain necessary "
                           "to achieve a new grain size. "
                           "List must have one more entry than the Phase transition depths.");
        prm.declare_entry ("Recrystallized grain size", "",
                           Patterns::List (Patterns::Double (0.)),
                           "The grain size $d_{ph}$ to that a phase will be reduced to when crossing a phase transition. "
                           "When set to zero, grain size will not be reduced. "
                           "List must have the same number of entries as Phase transition depths. "
                           "Units: \\si{\\meter}.");
        prm.declare_entry ("Phase volume fraction", "0.4",
                           Patterns::Double (0., 1.),
                           "The volume fraction of one of the phases in the two-phase damage model of Bercovici and Ricard (2012). "
                           "The volume fraction of the other phase can be simply calculated by subtracting from one. "
                           "This parameter is only used in the pinned state grain damage formulation."
                           "Units: none.");
        prm.declare_entry ("Grain size evolution formulation", "paleowattmeter",
                           Patterns::Selection ("paleowattmeter|paleopiezometer|pinned grain damage"),
                           "A flag indicating whether the material model should use the "
                           "paleowattmeter approach of Austin and Evans (2007) for grain size reduction "
                           "in the dislocation creep regime, the paleopiezometer approach "
                           "from Hall and Parmetier (2003), or the pinned grain damage approach "
                           "from Mulyukova and Bercovici (2018).");
        prm.declare_entry ("Use paleowattmeter", "default",
                           Patterns::Selection ("true|false|default"),
                           "A flag indicating whether the computation should use the "
                           "paleowattmeter approach of Austin and Evans (2007) for grain size reduction "
                           "in the dislocation creep regime (if true) or the paleopiezometer approach "
                           "from Hall and Parmetier (2003) (if false). This parameter has been removed. "
                           "Use 'Grain size evolution formulation' instead.");
        prm.declare_entry ("Average specific grain boundary energy", "1.0",
                           Patterns::List (Patterns::Double (0.)),
                           "The average specific grain boundary energy $\\gamma$. "
                           "List must have one more entry than the Phase transition depths. "
                           "Units: \\si{\\joule\\per\\meter\\squared}.");
        prm.declare_entry ("Work fraction for boundary area change", "0.1",
                           Patterns::List (Patterns::Double (0.)),
                           "The fraction $\\chi$ of work done by dislocation creep to change the grain boundary area. "
                           "List must have one more entry than the Phase transition depths. "
                           "Units: \\si{\\joule\\per\\meter\\squared}.");
        prm.declare_entry ("Geometric constant", "3.",
                           Patterns::List (Patterns::Double (0.)),
                           "The geometric constant $c$ used in the paleowattmeter grain size reduction law. "
                           "List must have one more entry than the Phase transition depths. "
                           "Units: none.");
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

        prm.enter_subsection("Grain damage partitioning");
        {
          prm.declare_entry ("Temperature for minimum grain damage partitioning", "1600",
                             Patterns::Double (0.),
                             "This parameter determines the temperature at which the computed coefficient of shear energy "
                             "partitioned into grain damage is minimum. This is used in the pinned state limit of the grain "
                             "size evolution. One choice of this parameter is the mantle temperature at the ridge axis, "
                             "see Mulyukova and Bercovici (2018) for details.");
          prm.declare_entry ("Temperature for maximum grain damage partitioning", "283",
                             Patterns::Double (0.),
                             "This parameter determines the temperature at which the computed coefficient of shear energy "
                             "partitioned into grain damage is maximum. This is used in the pinned state limit of the grain "
                             "size evolution. One choice of this parameter is the surface temperature of the seafloor, see "
                             "Mulyukova and Bercovici (2018) for details.");
          prm.declare_entry ("Minimum grain size reduction work fraction", "1e-12",
                             Patterns::Double (0., 1.),
                             "This parameter determines the minimum value of the partitioning coefficient, which governs "
                             "the amount of shear heating partitioned into grain damage in the pinned state limit.");
          prm.declare_entry ("Maximum grain size reduction work fraction", "1e-1",
                             Patterns::Double (0., 1.),
                             "This parameter determines the maximum value of the partitioning coefficient, which governs "
                             "the amount of shear heating partitioned into grain damage in the pinned state limit.");
          prm.declare_entry ("Grain size reduction work fraction exponent", "10",
                             Patterns::Double (0.),
                             "This parameter determines the variability in how much shear heating is partitioned into "
                             "grain damage. A higher value suggests a wider temperature range over which the partitioning "
                             "coefficient is high.");
        }
        prm.leave_subsection();
      }



      template <int dim>
      void
      GrainSizeEvolution<dim>::parse_parameters (ParameterHandler &prm)
      {
        AssertThrow (this->introspection().get_number_of_fields_of_type(CompositionalFieldDescription::grain_size) == 1,
                     ExcMessage("The 'grain size' material model only works if exactly one compositional "
                                "field with type 'grain size' is present. It looks like there are " +
                                std::to_string(this->introspection().get_number_of_fields_of_type(CompositionalFieldDescription::grain_size))
                                + " fields of this type."));

        recrystallized_grain_size = Utilities::string_to_double
                                    (Utilities::split_string_list(prm.get ("Recrystallized grain size")));

        if (recrystallized_grain_size.size() != n_phase_transitions)
          AssertThrow(false,
                      ExcMessage("Error: The list of recrystallized grain sizes has to have as many entries as there are phases."));

        for (unsigned int i=1; i<n_phase_transitions; ++i)
          AssertThrow(phase_function->get_transition_depth(i-1) < phase_function->get_transition_depth(i),
                      ExcMessage("Error: Phase transition depths have to be sorted in ascending order!"));

        // grain evolution parameters
        grain_growth_activation_energy        = Utilities::string_to_double
                                                (Utilities::split_string_list(prm.get ("Grain growth activation energy")));
        grain_growth_activation_volume        = Utilities::string_to_double
                                                (Utilities::split_string_list(prm.get ("Grain growth activation volume")));
        grain_growth_rate_constant            = Utilities::string_to_double
                                                (Utilities::split_string_list(prm.get ("Grain growth rate constant")));
        grain_growth_exponent                 = Utilities::string_to_double
                                                (Utilities::split_string_list(prm.get ("Grain growth exponent")));
        minimum_grain_size                    = prm.get_double("Minimum grain size");
        reciprocal_required_strain            = Utilities::string_to_double
                                                (Utilities::split_string_list(prm.get ("Reciprocal required strain")));

        grain_size_evolution_formulation      = Formulation::parse(prm.get("Grain size evolution formulation"));

        AssertThrow((grain_size_evolution_formulation != Formulation::paleopiezometer || !this->get_heating_model_manager().shear_heating_enabled()),
                    ExcMessage("Shear heating output should not be used with the Paleopiezometer grain damage formulation."));

        // TODO: Remove deprecated parameter in next release.
        const std::string use_paleowattmeter  = prm.get ("Use paleowattmeter");
        Assert(use_paleowattmeter == "default",
               ExcMessage("The parameter 'Use paleowattmeter' has been removed. "
                          "Use the parameter 'Grain size evolution formulation instead'."));

        const double volume_fraction_phase_one = prm.get_double ("Phase volume fraction");

        AssertThrow(volume_fraction_phase_one != 0. && volume_fraction_phase_one != 1.,
                    ExcMessage("Volume fraction must be between (0, 1) to use two phase damage in the pinned state!"));

        phase_distribution = phase_distribution_function(volume_fraction_phase_one);
        roughness_to_grain_size = roughness_to_grain_size_factor(volume_fraction_phase_one);

        grain_boundary_energy                 = Utilities::string_to_double
                                                (Utilities::split_string_list(prm.get ("Average specific grain boundary energy")));
        boundary_area_change_work_fraction    = Utilities::string_to_double
                                                (Utilities::split_string_list(prm.get ("Work fraction for boundary area change")));
        geometric_constant                    = Utilities::string_to_double
                                                (Utilities::split_string_list(prm.get ("Geometric constant")));

        if (grain_size_evolution_formulation == Formulation::pinned_grain_damage)
          {
            prm.enter_subsection("Grain damage partitioning");
            {
              grain_size_reduction_work_fraction_exponent = prm.get_double ("Grain size reduction work fraction exponent");
              maximum_grain_size_reduction_work_fraction  = prm.get_double ("Maximum grain size reduction work fraction");
              minimum_grain_size_reduction_work_fraction  = prm.get_double ("Minimum grain size reduction work fraction");

              AssertThrow(maximum_grain_size_reduction_work_fraction > 0. && maximum_grain_size_reduction_work_fraction < 1.,
                          ExcMessage("Maximum grain size reduction work fraction cannot be smaller or equal to 0 or larger or equal to 1."));
              AssertThrow(minimum_grain_size_reduction_work_fraction > 0. && minimum_grain_size_reduction_work_fraction < 1.,
                          ExcMessage("Minimum grain size reduction work fraction cannot be smaller or equal to 0 or larger or equal to 1."));
              AssertThrow(maximum_grain_size_reduction_work_fraction >= minimum_grain_size_reduction_work_fraction,
                          ExcMessage("Maximum grain size reduction work fraction must be larger than minimum grain size reduction work fraction."));

              const double temperature_minimum_partition  = prm.get_double ("Temperature for minimum grain damage partitioning");
              const double temperature_maximum_partition  = prm.get_double ("Temperature for maximum grain damage partitioning");

              AssertThrow(temperature_minimum_partition > temperature_maximum_partition,
                          ExcMessage("Temperature for minimum grain damage partitioning must be larger than Temperature for maximum grain damage partitioning."));

              temperature_minimum_partitioning_power = std::pow(temperature_minimum_partition,grain_size_reduction_work_fraction_exponent);
              temperature_maximum_partitioning_power = std::pow(temperature_maximum_partition,grain_size_reduction_work_fraction_exponent);
            }
            prm.leave_subsection();
          }

        // TODO: Remove deprecated parameters in next release.
        const double pv_grain_size_scaling         = prm.get_double ("Lower mantle grain size scaling");
        AssertThrow(pv_grain_size_scaling == 1.0,
                    ExcMessage("Error: The 'Lower mantle grain size scaling' parameter "
                               "has been removed. Please remove it from your input file. For models "
                               "with large spatial variations in grain size, please advect your "
                               "grain size on particles."));

        // TODO: Remove deprecated parameters in next release.
        const bool advect_log_grainsize            = prm.get_bool ("Advect logarithm of grain size");
        AssertThrow(advect_log_grainsize == false,
                    ExcMessage("Error: The 'Advect logarithm of grain size' parameter "
                               "has been removed. Please remove it from your input file. For models "
                               "with large spatial variations in grain size, please advect your "
                               "grain size on particles."));

        if (grain_growth_activation_energy.size() != grain_growth_activation_volume.size() ||
            grain_growth_activation_energy.size() != grain_growth_rate_constant.size() ||
            grain_growth_activation_energy.size() != grain_growth_exponent.size())
          AssertThrow(false,
                      ExcMessage("Error: The lists of grain size evolution and flow law parameters "
                                 "need to have the same length!"));

        if (grain_size_evolution_formulation == Formulation::paleowattmeter)
          {
            if (grain_growth_activation_energy.size() != grain_boundary_energy.size() ||
                grain_growth_activation_energy.size() != boundary_area_change_work_fraction.size() ||
                grain_growth_activation_energy.size() != geometric_constant.size() )
              AssertThrow(false,
                          ExcMessage("Error: One of the lists of grain size evolution parameters "
                                     "given for the paleowattmeter does not have the correct length!"));
          }
        else if (grain_size_evolution_formulation == Formulation::paleopiezometer)
          {
            AssertThrow(grain_growth_activation_energy.size() == reciprocal_required_strain.size(),
                        ExcMessage("Error: The list of grain size evolution parameters in the "
                                   "paleopiezometer does not have the correct length!"));
          }
        else if (grain_size_evolution_formulation == Formulation::pinned_grain_damage)
          {
            AssertThrow(n_phase_transitions == 0,
                        ExcMessage("Error: Currently, the pinned grain damage formulation is only implemented for one mineral phase."));
          }
        else
          AssertThrow(false,
                      ExcMessage("Error: The size of lists in grain size evolution and flow law parameters "
                                 "should follow either of the 'paleowattmeter|paleopiezometer|pinned grain damage' "
                                 "formulations!"));

        AssertThrow(grain_growth_activation_energy.size() == n_phase_transitions+1,
                    ExcMessage("Error: The lists of grain size evolution and flow law parameters need to "
                               "have exactly one more entry than the number of phase transitions "
                               "(which is defined by the length of the lists of phase transition depths, ...)!"));
      }



      template <int dim>
      void
      GrainSizeEvolution<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        // These properties will be used by the heating model to reduce
        // shear heating by the amount of work done to reduce grain size.
        if (out.template get_additional_output<HeatingModel::ShearHeatingOutputs<dim>>() == nullptr)
          {
            const unsigned int n_points = out.n_evaluation_points();
            out.additional_outputs.push_back(
              std::make_unique<HeatingModel::ShearHeatingOutputs<dim>> (n_points));
          }
      }



      template <int dim>
      void
      GrainSizeEvolution<dim>::fill_additional_outputs (const typename MaterialModel::MaterialModelInputs<dim> &in,
                                                        const typename MaterialModel::MaterialModelOutputs<dim> &out,
                                                        const std::vector<unsigned int> &phase_indices,
                                                        const std::vector<double> &dislocation_viscosities,
                                                        std::vector<std::unique_ptr<MaterialModel::AdditionalMaterialOutputs<dim>>> &additional_outputs) const
      {
        for (auto &additional_output: additional_outputs)
          if (HeatingModel::ShearHeatingOutputs<dim> *shear_heating_out = dynamic_cast<HeatingModel::ShearHeatingOutputs<dim> *>(additional_output.get()))
            {
              for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
                {
                  if (grain_size_evolution_formulation == Formulation::paleowattmeter)
                    {
                      const double f = boundary_area_change_work_fraction[phase_indices[i]];
                      shear_heating_out->shear_heating_work_fractions[i] = 1. - f * out.viscosities[i] / dislocation_viscosities[i];
                    }
                  else if (grain_size_evolution_formulation == Formulation::pinned_grain_damage)
                    {
                      const double f = compute_partitioning_fraction(in.temperature[i]);
                      shear_heating_out->shear_heating_work_fractions[i] = 1. - f;
                    }
                  else
                    AssertThrow(false, ExcNotImplemented());
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
  namespace ReactionModel \
  { \
    template class GrainSizeEvolution<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
