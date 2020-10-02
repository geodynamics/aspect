/*
  Copyright (C) 2019 - 2020 by the authors of the ASPECT code.

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


#include <aspect/material_model/rheology/strain_dependent.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>
#include <aspect/utilities.h>
#include <aspect/postprocess/particles.h>
#include <aspect/particle/property/interface.h>
#include <aspect/simulator.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      template <int dim>
      void
      StrainDependent<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Strain weakening mechanism", "default",
                           Patterns::Selection("none|finite strain tensor|total strain|plastic weakening with plastic strain only|plastic weakening with total strain only|plastic weakening with plastic strain and viscous weakening with viscous strain|viscous weakening with viscous strain only|default"),
                           "Whether to apply strain weakening to viscosity, cohesion and internal angle"
                           "of friction based on accumulated finite strain, and if yes, which method to "
                           "use. The following methods are available:"
                           "\n\n"
                           "\\item ``none'': No strain weakening is applied. "
                           "\n\n"
                           "\\item ``finite strain tensor'': The full finite strain tensor is tracked, "
                           "and its second invariant is used to weaken both the plastic yield stress "
                           "(specifically, the cohesion and friction angle) and the pre-yield viscosity "
                           "that arises from diffusion and/or dislocation creep."
                           "\n\n"
                           "\\item ``total strain'': The finite strain is approximated as the product of "
                           "the second invariant of the strain rate in each time step and the time step "
                           "size, and this quantity is integrated and tracked over time. It is used to "
                           "weaken both the plastic yield stress (specifically, the cohesion and friction "
                           "angle) and the pre-yield viscosity."
                           "\n\n"
                           "\\item ``plastic weakening with plastic strain only'': The finite strain is "
                           "approximated as the product of the second invariant of the strain rate"
                           "in each time step and the time step size in regions where material is "
                           "plastically yielding. This quantity is integrated and tracked over time, and "
                           "used to weaken the cohesion and friction angle. The pre-yield viscosity is "
                           "not weakened."
                           "\n\n"
                           "\\item ``plastic weakening with total strain only'': The finite strain is "
                           "approximated as the product of the second invariant of the strain rate in each "
                           "time step and the time step size, and this quantity is integrated and tracked "
                           "over time. It is used to weaken the plastic yield stress (specifically, the "
                           "cohesion and internal friction angle). The pre-yield viscosity is not weakened."
                           "\n\n"
                           "\\item ``plastic weakening with plastic strain and viscous weakening with viscous strain'': "
                           "Both the finite strain accumulated by plastic deformation and by viscous deformation are "
                           "computed separately (each approximated as the product of the second invariant of the "
                           "corresponding strain rate in each time step and the time step size). The plastic strain "
                           "is used to weaken the plastic yield stress (specifically, the cohesion and yield angle), and "
                           "the viscous strain is used to weaken the pre-yield viscosity."
                           "\n\n"
                           "\\item ``viscous weakening with viscous strain only'': The finite strain is "
                           "approximated as the product of the second invariant of the strain rate "
                           "in each time step and the time step size in regions where material is "
                           "not plastically yielding. This quantity is integrated and tracked over time, and "
                           "used to weaken the the pre-yield viscosity. The cohesion and friction angle are "
                           "not weakened."
                           "\n\n"
                           "\\item ``default'': The default option has the same behavior as ``none'', "
                           "but is there to make sure that the original parameters for specifying the "
                           "strain weakening mechanism (``Use plastic/viscous strain weakening'') are still allowed, "
                           "but to guarantee that one uses either the old parameter names or the new ones, "
                           "never both."
                           "\n\n"
                           "If a compositional field named 'noninitial\\_plastic\\_strain' is "
                           "included in the parameter file, this field will automatically be excluded from "
                           "from volume fraction calculation and track the cumulative plastic strain with "
                           "the initial plastic strain values removed.");

        prm.declare_entry ("Start plasticity strain weakening intervals", "0.",
                           Patterns::List(Patterns::Double (0.)),
                           "List of strain weakening interval initial strains "
                           "for the cohesion and friction angle parameters of the "
                           "background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value. Units: None.");

        prm.declare_entry ("End plasticity strain weakening intervals", "1.",
                           Patterns::List(Patterns::Double (0.)),
                           "List of strain weakening interval final strains "
                           "for the cohesion and friction angle parameters of the "
                           "background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value.  Units: None.");

        prm.declare_entry ("Cohesion strain weakening factors", "1.",
                           Patterns::List(Patterns::Double (0.)),
                           "List of cohesion strain weakening factors "
                           "for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value.  Units: None.");

        prm.declare_entry ("Friction strain weakening factors", "1.",
                           Patterns::List(Patterns::Double (0.)),
                           "List of friction strain weakening factors "
                           "for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value.  Units: None.");

        prm.declare_entry ("Start prefactor strain weakening intervals", "0.",
                           Patterns::List(Patterns::Double (0.)),
                           "List of strain weakening interval initial strains "
                           "for the diffusion and dislocation prefactor parameters of the "
                           "background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value.  Units: None.");

        prm.declare_entry ("End prefactor strain weakening intervals", "1.",
                           Patterns::List(Patterns::Double (0.)),
                           "List of strain weakening interval final strains "
                           "for the diffusion and dislocation prefactor parameters of the "
                           "background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value.  Units: None.");

        prm.declare_entry ("Prefactor strain weakening factors", "1.",
                           Patterns::List(Patterns::Double(0., 1.)),
                           "List of viscous strain weakening factors "
                           "for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "If only one value is given, then all use the same value.  Units: None.");

        prm.declare_entry ("Strain healing mechanism", "no healing",
                           Patterns::Selection("no healing|temperature dependent"),
                           "Whether to apply strain healing to plastic yielding and viscosity terms, "
                           "and if yes, which method to use. The following methods are available:"
                           "\n\n"
                           "\\item ``no healing'': No strain healing is applied. "
                           "\n\n"
                           "\\item ``temperature dependent'': Purely temperature dependent "
                           "strain healing applied to plastic yielding and viscosity terms, similar "
                           "to the temperature-dependent Frank Kamenetskii formulation, computes "
                           "strain healing as removing strain as a function of temperature, time, "
                           "and a user-defined healing rate and prefactor "
                           "as done in Fuchs and Becker, 2019, for mantle convection");

        prm.declare_entry ("Strain healing temperature dependent recovery rate", "1.e-15", Patterns::Double(0),
                           "Recovery rate prefactor for temperature dependent "
                           "strain healing. Units: $1/s$");

        prm.declare_entry ("Strain healing temperature dependent prefactor", "15.", Patterns::Double(0),
                           "Prefactor for temperature dependent "
                           "strain healing. Units: None");
      }

      template <int dim>
      void
      StrainDependent<dim>::parse_parameters (ParameterHandler &prm)
      {
        // Get the number of fields for composition-dependent material properties
        const unsigned int n_fields = this->n_compositional_fields() + 1;

        // number of required compositional fields for full finite strain tensor
        const unsigned int s = Tensor<2,dim>::n_independent_components;

        // Strain weakening parameters
        if (prm.get ("Strain weakening mechanism") == "none")
          weakening_mechanism = none;
        else if (prm.get ("Strain weakening mechanism") == "finite strain tensor")
          weakening_mechanism = finite_strain_tensor;
        else if (prm.get ("Strain weakening mechanism") == "total strain")
          weakening_mechanism = total_strain;
        else if (prm.get ("Strain weakening mechanism") == "plastic weakening with plastic strain only")
          weakening_mechanism = plastic_weakening_with_plastic_strain_only;
        else if (prm.get ("Strain weakening mechanism") == "plastic weakening with total strain only")
          weakening_mechanism = plastic_weakening_with_total_strain_only;
        else if (prm.get ("Strain weakening mechanism") == "plastic weakening with plastic strain and viscous weakening with viscous strain")
          weakening_mechanism = plastic_weakening_with_plastic_strain_and_viscous_weakening_with_viscous_strain;
        else if (prm.get ("Strain weakening mechanism") == "viscous weakening with viscous strain only")
          weakening_mechanism = viscous_weakening_with_viscous_strain_only;
        else if (prm.get ("Strain weakening mechanism") == "default")
          weakening_mechanism = none;
        else
          AssertThrow(false, ExcMessage("Not a valid Strain weakening mechanism!"));

        if (weakening_mechanism == plastic_weakening_with_plastic_strain_only
            || weakening_mechanism == plastic_weakening_with_plastic_strain_and_viscous_weakening_with_viscous_strain)
          {
            AssertThrow(this->introspection().compositional_name_exists("plastic_strain"),
                        ExcMessage("Material model visco_plastic with plastic strain weakening only works if there is a "
                                   "compositional field called plastic_strain."));
          }

        if (weakening_mechanism == viscous_weakening_with_viscous_strain_only
            || weakening_mechanism == plastic_weakening_with_plastic_strain_and_viscous_weakening_with_viscous_strain)
          {
            AssertThrow(this->introspection().compositional_name_exists("viscous_strain"),
                        ExcMessage("Material model visco_plastic with viscous strain weakening only works if there is a "
                                   "compositional field called viscous_strain."));
          }

        if (weakening_mechanism == finite_strain_tensor)
          {
            AssertThrow(this->n_compositional_fields() >= s,
                        ExcMessage("There must be enough compositional fields to track all components of the finite strain tensor (4 in 2D, 9 in 3D). "));
            // Assert that fields exist and that they are in the right order
            const unsigned int n_s11 = this->introspection().compositional_index_for_name("s11");
            const unsigned int n_s12 = this->introspection().compositional_index_for_name("s12");
            const unsigned int n_s21 = this->introspection().compositional_index_for_name("s21");
            const unsigned int n_s22 = this->introspection().compositional_index_for_name("s22");
            if (dim==2)
              {
                AssertThrow(n_s11 < n_s12 && n_s12 < n_s21 && n_s21 < n_s22,
                            ExcMessage("A material model with strain weakening using the full strain tensor only works if there are "
                                       "compositional fields called sij, with i=1,..,dim and j=1,...,dim listed in the following order: "
                                       "s11, s12, s21, s22."));
                AssertThrow(n_s22 == n_s11+s-1, ExcMessage("The strain tensor components should be represented by consecutive fields."))
              }
            if (dim==3)
              {
                const unsigned int n_s13 = this->introspection().compositional_index_for_name("s13");
                const unsigned int n_s23 = this->introspection().compositional_index_for_name("s23");
                const unsigned int n_s31 = this->introspection().compositional_index_for_name("s31");
                const unsigned int n_s32 = this->introspection().compositional_index_for_name("s32");
                const unsigned int n_s33 = this->introspection().compositional_index_for_name("s33");
                AssertThrow(n_s11 < n_s12 && n_s12 < n_s13 && n_s13 < n_s21 && n_s21 < n_s22 && n_s22 < n_s23 && n_s23 < n_s31 && n_s31 < n_s32 && n_s32 < n_s33,
                            ExcMessage("A material model with strain weakening using the full strain tensor only works if there are "
                                       "compositional fields called sij, with i=1,..,dim and j=1,...,dim listed in the following order: "
                                       "s11, s12, s13, s21, s22, s23, s31, s32, s33."));
                AssertThrow(n_s33 == n_s11+s-1, ExcMessage("The strain tensor components should be represented by consecutive fields."));
              }
          }


        if (weakening_mechanism == total_strain || weakening_mechanism == plastic_weakening_with_total_strain_only)
          AssertThrow(this->introspection().compositional_name_exists("total_strain"),
                      ExcMessage("Material model visco_plastic with total strain weakening only works if there is a "
                                 "compositional field called total_strain."));

        // Currently, it only makes sense to use this material model with strain weakening when the
        // nonlinear solver scheme does a single advection iteration. More than one nonlinear advection
        // iteration will result in the incorrect value of strain being used in the material model, as
        // the compositional fields representing strain are updated through the reaction terms.
        if (weakening_mechanism != none)
          {
            AssertThrow((this->get_parameters().nonlinear_solver ==
                         Parameters<dim>::NonlinearSolver::single_Advection_single_Stokes
                         ||
                         this->get_parameters().nonlinear_solver ==
                         Parameters<dim>::NonlinearSolver::single_Advection_iterated_Stokes
                         ||
                         this->get_parameters().nonlinear_solver ==
                         Parameters<dim>::NonlinearSolver::single_Advection_iterated_Newton_Stokes),
                        ExcMessage("The material model will only work with the nonlinear "
                                   "solver schemes 'single Advection, single Stokes', "
                                   "'single Advection, iterated Stokes', and "
                                   "'single Advection, iterated_Newton_Stokes' when strain "
                                   "weakening is enabled."));
          }



        start_plastic_strain_weakening_intervals = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Start plasticity strain weakening intervals"))),
                                                   n_fields,
                                                   "Start plasticity strain weakening intervals");

        end_plastic_strain_weakening_intervals = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("End plasticity strain weakening intervals"))),
                                                 n_fields,
                                                 "End plasticity strain weakening intervals");

        start_viscous_strain_weakening_intervals = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Start prefactor strain weakening intervals"))),
                                                   n_fields,
                                                   "Start prefactor strain weakening intervals");

        end_viscous_strain_weakening_intervals = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("End prefactor strain weakening intervals"))),
                                                 n_fields,
                                                 "End prefactor strain weakening intervals");

        viscous_strain_weakening_factors = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Prefactor strain weakening factors"))),
                                                                                   n_fields,
                                                                                   "Prefactor strain weakening factors");

        cohesion_strain_weakening_factors = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Cohesion strain weakening factors"))),
                                                                                    n_fields,
                                                                                    "Cohesion strain weakening factors");

        friction_strain_weakening_factors = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Friction strain weakening factors"))),
                                                                                    n_fields,
                                                                                    "Friction strain weakening factors");

        if (prm.get ("Strain healing mechanism") == "no healing")
          healing_mechanism = no_healing;
        else if (prm.get ("Strain healing mechanism") == "temperature dependent")
          healing_mechanism = temperature_dependent;
        else
          AssertThrow(false, ExcMessage("Not a valid Strain healing mechanism!"));

        // Currently this functionality only works in field composition
        if (healing_mechanism != no_healing && this->get_postprocess_manager().template has_matching_postprocessor<Postprocess::Particles<dim> >())
          {
            const Postprocess::Particles<dim> &particle_postprocessor = this->get_postprocess_manager().template get_matching_postprocessor<Postprocess::Particles<dim> >();
            const Particle::Property::Manager<dim> &particle_property_manager = particle_postprocessor.get_particle_world().get_property_manager();

            AssertThrow(particle_property_manager.plugin_name_exists("viscoplastic strain invariants") == false, ExcMessage("This healing mechanism currently does not work if the strain is tracked on particles."));
          }

        // Temperature dependent strain healing requires that adiabatic surface temperature is non zero
        if (healing_mechanism == temperature_dependent)
          {
            AssertThrow (this->get_adiabatic_surface_temperature() > 0.0,
                         ExcMessage("The temperature dependent strain healing can only be used when the adiabatic "
                                    "surface temperature (reference_temperature in equation for strain healing) "
                                    "is non-zero."));
          }

        strain_healing_temperature_dependent_recovery_rate = prm.get_double ("Strain healing temperature dependent recovery rate");

        strain_healing_temperature_dependent_prefactor = prm.get_double ("Strain healing temperature dependent prefactor");
      }


      template <int dim>
      std::array<double, 3>
      StrainDependent<dim>::
      compute_strain_weakening_factors(const unsigned int j,
                                       const std::vector<double> &composition) const
      {
        double viscous_weakening = 1.0;
        std::pair<double, double> brittle_weakening (1.0, 1.0);

        switch (weakening_mechanism)
          {
            case none:
            {
              break;
            }
            case finite_strain_tensor:
            {
              // Calculate second invariant of left stretching tensor "L"
              Tensor<2,dim> strain;
              for (unsigned int q = 0; q < Tensor<2,dim>::n_independent_components ; ++q)
                strain[Tensor<2,dim>::unrolled_to_component_indices(q)] = composition[q];
              const SymmetricTensor<2,dim> L = symmetrize( strain * transpose(strain) );

              const double strain_ii = std::fabs(second_invariant(L));
              brittle_weakening = calculate_plastic_weakening(strain_ii, j);
              viscous_weakening = calculate_viscous_weakening(strain_ii, j);
              break;
            }
            case total_strain:
            {
              const unsigned int total_strain_index = this->introspection().compositional_index_for_name("total_strain");
              brittle_weakening = calculate_plastic_weakening(composition[total_strain_index], j);
              viscous_weakening = calculate_viscous_weakening(composition[total_strain_index], j);
              break;
            }
            case plastic_weakening_with_total_strain_only:
            {
              const unsigned int total_strain_index = this->introspection().compositional_index_for_name("total_strain");
              brittle_weakening = calculate_plastic_weakening(composition[total_strain_index], j);
              break;
            }
            case plastic_weakening_with_plastic_strain_only:
            {
              const unsigned int plastic_strain_index = this->introspection().compositional_index_for_name("plastic_strain");
              brittle_weakening = calculate_plastic_weakening(composition[plastic_strain_index], j);
              break;
            }
            case plastic_weakening_with_plastic_strain_and_viscous_weakening_with_viscous_strain:
            {
              const unsigned int plastic_strain_index = this->introspection().compositional_index_for_name("plastic_strain");
              brittle_weakening = calculate_plastic_weakening(composition[plastic_strain_index], j);
              const unsigned int viscous_strain_index = this->introspection().compositional_index_for_name("viscous_strain");
              viscous_weakening = calculate_viscous_weakening(composition[viscous_strain_index], j);
              break;
            }
            case viscous_weakening_with_viscous_strain_only:
            {
              const unsigned int viscous_strain_index = this->introspection().compositional_index_for_name("viscous_strain");
              viscous_weakening = calculate_viscous_weakening(composition[viscous_strain_index], j);
              break;
            }
            default:
            {
              AssertThrow(false, ExcNotImplemented());
              break;
            }
          }

        const std::array<double, 3> weakening_factors = {{brittle_weakening.first,brittle_weakening.second,viscous_weakening}};

        return weakening_factors;

      }

      template <int dim>
      double
      StrainDependent<dim>::
      calculate_strain_healing(const MaterialModel::MaterialModelInputs<dim> &in,
                               const unsigned int j) const
      {
        const double reference_temperature = this->get_adiabatic_surface_temperature();
        double healed_strain = 0.0;

        switch (healing_mechanism)
          {
            case no_healing:
            {
              break;
            }
            case temperature_dependent:
            {
              healed_strain = strain_healing_temperature_dependent_recovery_rate *
                              std::exp(-strain_healing_temperature_dependent_prefactor * 0.5 * (1.0 - in.temperature[j]/reference_temperature))
                              * this->get_timestep();
              break;
            }
          }
        return healed_strain;
      }

      template <int dim>
      std::pair<double, double>
      StrainDependent<dim>::
      calculate_plastic_weakening(const double strain_ii,
                                  const unsigned int j) const
      {
        // Constrain the second strain invariant of the previous timestep by the strain interval
        const double cut_off_strain_ii = std::max(std::min(strain_ii,end_plastic_strain_weakening_intervals[j]),start_plastic_strain_weakening_intervals[j]);

        // Linear strain weakening of cohesion and internal friction angle between specified strain values
        const double strain_fraction = (cut_off_strain_ii - start_plastic_strain_weakening_intervals[j]) /
                                       (start_plastic_strain_weakening_intervals[j] - end_plastic_strain_weakening_intervals[j]);

        const double weakening_cohesion = 1. + (1. - cohesion_strain_weakening_factors[j]) * strain_fraction;
        const double weakening_friction = 1. + (1. - friction_strain_weakening_factors[j]) * strain_fraction;

        return std::make_pair (weakening_cohesion, weakening_friction);
      }


      template <int dim>
      double
      StrainDependent<dim>::
      calculate_viscous_weakening(const double strain_ii,
                                  const unsigned int j) const
      {
        // Constrain the second strain invariant of the previous timestep by the strain interval
        const double cut_off_strain_ii = std::max(std::min(strain_ii,end_viscous_strain_weakening_intervals[j]),start_viscous_strain_weakening_intervals[j]);

        // Linear strain weakening of the viscous flow law prefactors between specified strain values
        const double strain_fraction = (cut_off_strain_ii - start_viscous_strain_weakening_intervals[j]) /
                                       (start_viscous_strain_weakening_intervals[j] - end_viscous_strain_weakening_intervals[j]);
        return 1. + ( 1. - viscous_strain_weakening_factors[j] ) * strain_fraction;
      }

      template <int dim>
      void
      StrainDependent<dim>::
      fill_reaction_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
                             const int i,
                             const double min_strain_rate,
                             const bool plastic_yielding,
                             MaterialModel::MaterialModelOutputs<dim> &out) const
      {

        // If strain weakening is used, overwrite the first reaction term,
        // which represents the second invariant of the (plastic) strain tensor.
        // If plastic strain is tracked (so not the total strain), only overwrite
        // when plastically yielding.
        // If viscous strain is also tracked, overwrite the second reaction term as well.
        // Calculate changes in strain and update the reaction terms
        if  (this->simulator_is_past_initialization() && this->get_timestep_number() > 0 && in.requests_property(MaterialProperties::reaction_terms))
          {
            const double edot_ii = std::max(sqrt(std::fabs(second_invariant(deviator(in.strain_rate[i])))),min_strain_rate);
            double delta_e_ii = edot_ii*this->get_timestep();

            // Adjusting strain values to account for strain healing without exceeding an unreasonable range
            if (healing_mechanism != no_healing)
              {
                // Never heal more strain than exists
                delta_e_ii -= calculate_strain_healing(in,i);
              }
            if (weakening_mechanism == plastic_weakening_with_plastic_strain_only && plastic_yielding == true)
              out.reaction_terms[i][this->introspection().compositional_index_for_name("plastic_strain")] =
                std::max(delta_e_ii, -in.composition[i][this->introspection().compositional_index_for_name("plastic_strain")]);
            if (weakening_mechanism == viscous_weakening_with_viscous_strain_only && plastic_yielding == false)
              out.reaction_terms[i][this->introspection().compositional_index_for_name("viscous_strain")] =
                std::max(delta_e_ii, -in.composition[i][this->introspection().compositional_index_for_name("viscous_strain")]);
            if (weakening_mechanism == total_strain || weakening_mechanism == plastic_weakening_with_total_strain_only)
              out.reaction_terms[i][this->introspection().compositional_index_for_name("total_strain")] =
                std::max(delta_e_ii, -in.composition[i][this->introspection().compositional_index_for_name("total_strain")]);
            if (weakening_mechanism == plastic_weakening_with_plastic_strain_and_viscous_weakening_with_viscous_strain)
              {
                if (plastic_yielding == true)
                  out.reaction_terms[i][this->introspection().compositional_index_for_name("plastic_strain")] =
                    std::max(delta_e_ii, -in.composition[i][this->introspection().compositional_index_for_name("plastic_strain")]);
                else
                  out.reaction_terms[i][this->introspection().compositional_index_for_name("viscous_strain")] =
                    std::max(delta_e_ii, -in.composition[i][this->introspection().compositional_index_for_name("viscous_strain")]);
              }
            if (this->introspection().compositional_name_exists("noninitial_plastic_strain") && plastic_yielding == true)
              out.reaction_terms[i][this->introspection().compositional_index_for_name("noninitial_plastic_strain")] =
                std::max(delta_e_ii, -in.composition[i][this->introspection().compositional_index_for_name("noninitial_plastic_strain")]);
          }
      }


      template <int dim>
      void
      StrainDependent<dim>::
      compute_finite_strain_reaction_terms(const MaterialModel::MaterialModelInputs<dim> &in,
                                           MaterialModel::MaterialModelOutputs<dim> &out) const
      {

        if (in.current_cell.state() == IteratorState::valid && this->get_timestep_number() > 0 && in.requests_property(MaterialProperties::reaction_terms))
          {
            // We need the velocity gradient for the finite strain (they are not
            // in material model inputs), so we get them from the finite element.
            std::vector<Point<dim> > quadrature_positions(in.n_evaluation_points());
            for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
              quadrature_positions[i] = this->get_mapping().transform_real_to_unit_cell(in.current_cell, in.position[i]);

            FEValues<dim> fe_values (this->get_mapping(),
                                     this->get_fe(),
                                     Quadrature<dim>(quadrature_positions),
                                     update_gradients);

            std::vector<Tensor<2,dim> > velocity_gradients (quadrature_positions.size(), Tensor<2,dim>());

            fe_values.reinit (in.current_cell);
            fe_values[this->introspection().extractors.velocities].get_function_gradients (this->get_solution(),
                                                                                           velocity_gradients);

            // Assign the strain components to the compositional fields reaction terms.
            // If there are too many fields, we simply fill only the first fields with the
            // existing strain tensor components.

            for (unsigned int q=0; q < in.n_evaluation_points(); ++q)
              {
                if (in.current_cell.state() == IteratorState::valid && weakening_mechanism == finite_strain_tensor
                    && this->get_timestep_number() > 0 && in.requests_property(MaterialProperties::reaction_terms))

                  {
                    // Convert the compositional fields into the tensor quantity they represent.
                    Tensor<2,dim> strain;
                    const unsigned int n_first = this->introspection().compositional_index_for_name("s11");
                    for (unsigned int i = n_first; i < n_first + Tensor<2,dim>::n_independent_components ; ++i)
                      {
                        strain[Tensor<2,dim>::unrolled_to_component_indices(i)] = in.composition[q][i];
                      }

                    // Compute the strain accumulated in this timestep.
                    const Tensor<2,dim> strain_increment = this->get_timestep() * (velocity_gradients[q] * strain);

                    // Output the strain increment component-wise to its respective compositional field's reaction terms.
                    for (unsigned int i = n_first; i < n_first + Tensor<2,dim>::n_independent_components ; ++i)
                      {
                        out.reaction_terms[q][i] = strain_increment[Tensor<2,dim>::unrolled_to_component_indices(i)];
                      }
                  }
              }
          }
      }

      template <int dim>
      ComponentMask
      StrainDependent<dim>::
      get_strain_composition_mask() const
      {

        // Store which components to exclude during volume fraction computation.
        ComponentMask strain_mask(this->n_compositional_fields(),true);

        if (weakening_mechanism != none)
          {
            if (weakening_mechanism == plastic_weakening_with_plastic_strain_only || weakening_mechanism == plastic_weakening_with_plastic_strain_and_viscous_weakening_with_viscous_strain)
              strain_mask.set(this->introspection().compositional_index_for_name("plastic_strain"),false);

            if (weakening_mechanism == viscous_weakening_with_viscous_strain_only || weakening_mechanism == plastic_weakening_with_plastic_strain_and_viscous_weakening_with_viscous_strain)
              strain_mask.set(this->introspection().compositional_index_for_name("viscous_strain"),false);

            if (weakening_mechanism == total_strain || weakening_mechanism == plastic_weakening_with_total_strain_only)
              strain_mask.set(this->introspection().compositional_index_for_name("total_strain"),false);

            if (weakening_mechanism == finite_strain_tensor)
              {
                const unsigned int n_start = this->introspection().compositional_index_for_name("s11");
                for (unsigned int i = n_start; i < n_start + Tensor<2,dim>::n_independent_components ; ++i)
                  strain_mask.set(i,false);
              }
          }

        if (this->introspection().compositional_name_exists("noninitial_plastic_strain"))
          strain_mask.set(this->introspection().compositional_index_for_name("noninitial_plastic_strain"),false);

        return strain_mask;
      }

      template <int dim>
      WeakeningMechanism
      StrainDependent<dim>::
      get_weakening_mechanism() const
      {
        return weakening_mechanism;
      }

      template <int dim>
      HealingMechanism
      StrainDependent<dim>::
      get_healing_mechanism() const
      {
        return healing_mechanism;
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
    template class StrainDependent<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}

