/*
  Copyright (C) 2023 - 2024 by the authors of the ASPECT code.

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

#include <aspect/material_model/reactive_fluid_transport.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/numerics/fe_field_function.h>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    bool
    ReactiveFluidTransport<dim>::
    is_compressible () const
    {
      return base_model->is_compressible();
    }



    template <int dim>
    double
    ReactiveFluidTransport<dim>::
    reference_darcy_coefficient () const
    {
      if (fluid_solid_reaction_scheme == katz2003)
        {
          return katz2003_model.reference_darcy_coefficient();
        }
      else
        {
          // 0.01 = 1% melt
          return reference_permeability * Utilities::fixed_power<3>(0.01) / eta_f;
        }
    }



    template <int dim>
    double
    ReactiveFluidTransport<dim>::
    compute_bulk_density (const double porosity,
                          const double solid_density,
                          const double fluid_density) const
    {
      return (1 - porosity) * solid_density + porosity * fluid_density;
    }



    template <int dim>
    double
    ReactiveFluidTransport<dim>::
    compute_mass_fraction (const double volume_frac,
                           const double material_density,
                           const double bulk_density) const
    {
      return volume_frac * material_density / bulk_density;
    }



    template <int dim>
    void
    ReactiveFluidTransport<dim>::
    melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                    std::vector<double> &melt_fractions,
                    const MaterialModel::MaterialModelOutputs<dim> *out) const
    {
      const std::shared_ptr<const MeltOutputs<dim>> fluid_out = out->template get_additional_output_object<MeltOutputs<dim>>();
      if (fluid_out != nullptr && in.requests_property(MaterialProperties::additional_outputs))
        {
          for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
            {
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
              switch (fluid_solid_reaction_scheme)
                {
                  case no_reaction:
                  {
                    Assert(out != nullptr,
                           ExcMessage("The MaterialModelOutputs provided to the No Reactions reaction model is a nullptr."));
                    const double volume_fraction_porosity = in.composition[q][porosity_idx];
                    // Melt fractions outputs the mass fraction of equilibrium free water at each quadrature point. However,
                    // the free water is a volume fraction, while the bound water is a mass fraction. Convert the free water
                    // to a mass fraction so that we correctly compare the two values.
                    const std::shared_ptr<const MeltOutputs<dim>> fluid_out = out->template get_additional_output_object<MeltOutputs<dim>>();
                    const double bulk_density = compute_bulk_density(volume_fraction_porosity, out->densities[q], fluid_out->fluid_densities[q]);
                    const double mass_fraction_porosity = compute_mass_fraction(volume_fraction_porosity, fluid_out->fluid_densities[q], bulk_density);
                    // No reactions occur between the solid and fluid phases,
                    // and the fluid mass fraction (stored in the melt_fractions
                    // vector) is equal to the mass fraction of the porosity.
                    melt_fractions[q] = mass_fraction_porosity;
                    break;
                  }
                  case zero_solubility:
                  {
                    // The fluid volume fraction in equilibrium with the solid
                    // at any point (stored in the melt_fractions vector) is
                    // equal to the sum of the bound fluid content and porosity.
                    Assert(out != nullptr,
                           ExcMessage("The MaterialModelOutputs provided to the Zero Solubility reaction model is a nullptr."));
                    const double volume_fraction_porosity = in.composition[q][porosity_idx];
                    const unsigned int bound_fluid_idx = this->introspection().compositional_index_for_name("bound_fluid");
                    // Melt fractions outputs the mass fraction of equilibrium free water at each quadrature point. However,
                    // the free water is a volume fraction, while the bound water is a mass fraction. Convert the free water
                    // to a mass fraction so that we correctly compare the two values.
                    const std::shared_ptr<const MeltOutputs<dim>> fluid_out = out->template get_additional_output_object<MeltOutputs<dim>>();
                    const double bulk_density = compute_bulk_density(volume_fraction_porosity, out->densities[q], fluid_out->fluid_densities[q]);
                    const double mass_fraction_porosity = compute_mass_fraction(volume_fraction_porosity, fluid_out->fluid_densities[q], bulk_density);
                    melt_fractions[q] = in.composition[q][bound_fluid_idx] + mass_fraction_porosity;
                    break;
                  }
                  case tian_approximation:
                  {
                    Assert(out != nullptr,
                           ExcMessage("The MaterialModelOutputs provided to the Tian 2019 reaction model is a nullptr."));
                    // The tian2019_model melt_fraction function requires the mass fraction of free water at each quadrature
                    // point as an input. However, the free water is stored in the porosity field as a volume fraction.
                    // We therefore convert the free water to a mass fraction so that we have the correct input.
                    const double volume_fraction_porosity = in.composition[q][porosity_idx];
                    const std::shared_ptr<const MeltOutputs<dim>> fluid_out = out->template get_additional_output_object<MeltOutputs<dim>>();
                    const double bulk_density = compute_bulk_density(volume_fraction_porosity, out->densities[q], fluid_out->fluid_densities[q]);;
                    const double mass_fraction_porosity = compute_mass_fraction(volume_fraction_porosity, fluid_out->fluid_densities[q], bulk_density);
                    melt_fractions[q] = tian2019_model.melt_fraction(in, mass_fraction_porosity, q);
                    break;
                  }
                  case katz2003:
                  {
                    melt_fractions[q] = katz2003_model.melt_fraction(in.temperature[q],
                                                                     this->get_adiabatic_conditions().pressure(in.position[q]));
                    break;
                  }
                  default:
                  {
                    AssertThrow(false, ExcNotImplemented());
                    break;
                  }
                }
            }
        }
    }



    template <int dim>
    void
    ReactiveFluidTransport<dim>::initialize()
    {
      base_model->initialize();
    }



    template <int dim>
    void
    ReactiveFluidTransport<dim>::update()
    {
      base_model->update();
    }



    template <int dim>
    void
    ReactiveFluidTransport<dim>::evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                                          typename Interface<dim>::MaterialModelOutputs &out) const
    {
      base_model->evaluate(in,out);

      if (fluid_solid_reaction_scheme != katz2003)
        {
          const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

          // Modify the viscosity from the base model based on the presence of fluid.
          if (in.requests_property(MaterialProperties::viscosity))
            {
              // Scale the base model viscosity value based on the porosity.
              for (unsigned int q=0; q<out.n_evaluation_points(); ++q)
                {
                  const double porosity = std::max(in.composition[q][porosity_idx],0.0);
                  out.viscosities[q] *= (1.0 - porosity) * std::exp(- alpha_phi * porosity);
                }
            }

          // Fill the melt outputs if they exist. Note that the MeltOutputs class was originally
          // designed for two-phase flow material models in ASPECT that model the flow of melt,
          // but can be reused for a geofluid of arbitrary composition.
          const std::shared_ptr<MeltOutputs<dim>> fluid_out
            = out.template get_additional_output_object<MeltOutputs<dim>>();


          if (fluid_out != nullptr && in.requests_property(MaterialProperties::additional_outputs))
            {
              for (unsigned int q=0; q<out.n_evaluation_points(); ++q)
                {
                  double porosity = std::max(in.composition[q][porosity_idx],0.0);

                  fluid_out->fluid_viscosities[q] = eta_f;
                  fluid_out->permeabilities[q] = reference_permeability * Utilities::fixed_power<3>(porosity) * Utilities::fixed_power<2>(1.0-porosity);

                  fluid_out->fluid_densities[q] = reference_rho_f * std::exp(fluid_compressibility * (in.pressure[q] - this->get_surface_pressure()));
                  if (in.requests_property(MaterialProperties::viscosity))
                    {
                      const double phi_0 = 0.05;

                      // Limit the porosity to be no smaller than 1e-8 when
                      // calculating fluid effects on viscosities.
                      porosity = std::max(porosity,1e-8);
                      fluid_out->compaction_viscosities[q] = std::max(std::min(out.viscosities[q] * shear_to_bulk_viscosity_ratio * phi_0/porosity, max_compaction_visc), min_compaction_visc);
                    }
                }
            }

          const std::shared_ptr<ReactionRateOutputs<dim>> reaction_rate_out
            = out.template get_additional_output_object<ReactionRateOutputs<dim>>();

          if (this->get_parameters().use_operator_splitting && reaction_rate_out != nullptr
              && in.requests_property(MaterialProperties::reaction_rates))
            {
              // Fill reaction rate outputs if the model uses operator splitting.
              // Specifically, change the porosity (representing the amount of free fluid)
              // based on the water solubility and the fluid content.
              double reaction_time_step_size = 1.0;
              double reaction_fraction       = 0.0;
              if (this->simulator_is_past_initialization())
                {
                  const unsigned int number_of_reaction_steps = std::max(static_cast<unsigned int>(this->get_timestep() / this->get_parameters().reaction_time_step),
                                                                         std::max(this->get_parameters().reaction_steps_per_advection_step,1U));
                  reaction_time_step_size = this->get_timestep() / static_cast<double>(number_of_reaction_steps);
                  reaction_fraction       = reaction_time_step_size / fluid_reaction_time_scale;
                }

              std::vector<double> eq_free_fluid_fractions(out.n_evaluation_points());
              melt_fractions(in, eq_free_fluid_fractions, &out);

              for (unsigned int q=0; q<out.n_evaluation_points(); ++q)
                {
                  const unsigned int bound_fluid_idx = this->introspection().compositional_index_for_name("bound_fluid");

                  const double solid_density = out.densities[q];
                  const double fluid_density = fluid_out->fluid_densities[q];

                  const double volume_fraction_porosity = in.composition[q][porosity_idx];
                  const double bulk_density = compute_bulk_density(volume_fraction_porosity, solid_density, fluid_density);

                  // We need to convert the porosity, which is a volume fraction, to a mass fraction
                  // to determine the total mass fraction of water (across both phases). This is then
                  // used to ensure conservation of water mass during reactions.
                  const double mass_fraction_porosity = compute_mass_fraction(volume_fraction_porosity, fluid_density, bulk_density);

                  // The total mass fraction of water, which must be conserved after the reactions are done.
                  // This is given by:
                  // TM = Sm_wt * Sm_mass + Fl_wt * Fl_mass
                  // Where TM is the total mass of water, Sm_wt is the wt% of water in the solid, Sm_mass is the
                  // mass fraction of water in the solid, Fl_wt is the wt% of water in the fluid, and Fl_mass is
                  // the mass fraction of water in the fluid. Fl_wt is always 100%, because the fluid is assumed to
                  // be composed of only water. The bound_fluid composition gives Sm_wt, and we compute Sm_mass as
                  // was done when defining mass_frac_porosity.
                  const double total_mass_fraction_water = in.composition[q][bound_fluid_idx] * (1 - volume_fraction_porosity) * solid_density / bulk_density +
                                                           mass_fraction_porosity;

                  // Determine the change in the mass fraction of porosity from the equilibrium fluid fraction
                  // and update the mass fraction of porosity.
                  const double mass_fraction_porosity_change = (eq_free_fluid_fractions[q] - mass_fraction_porosity) * reaction_fraction;
                  const double new_mass_fraction_porosity = mass_fraction_porosity + mass_fraction_porosity_change;

                  // Since porosity is a volume fraction, convert the mass fraction change to a volume fraction change
                  // to update the porosity value. We cannot use the bulk density to do this, because we do not yet know
                  // the new bulk density of the rock after the reactions have been applied. We do this by substituting
                  // the  new bulk density as a function of the new volume fraction porosity into this equation:
                  // new_Fl_mass * new_bulk_density = new_volume_frac * fluid_density
                  // and solving for the new volume fraction. This gives us the following equation:
                  const double new_volume_fraction_porosity = solid_density * new_mass_fraction_porosity /
                                                              (fluid_density + solid_density * new_mass_fraction_porosity - fluid_density * new_mass_fraction_porosity);

                  // Now determine the new bulk density using the new volume fraction of porosity
                  const double new_bulk_density = compute_bulk_density(new_volume_fraction_porosity, solid_density, fluid_density);

                  // Determine the bound fluid mass fraction using the new volume fraction of porosity and the original
                  // mass fraction of water. We do this by determining the mass fraction of the solid phase:
                  const double mass_fraction_solid = compute_mass_fraction(1.0 - new_volume_fraction_porosity, solid_density, new_bulk_density);
                  // And now we use the total mass fraction of water determined above to ensure that we are conserving
                  // mass during the reactions. Sm_mass is then given by:
                  const double new_bound_fluid_mass_fraction = (total_mass_fraction_water - new_mass_fraction_porosity) /
                                                               (mass_fraction_solid);

                  // Determine the change in both the bound and free water as a mass fraction and volume fraction, respectively
                  double bound_fluid_change = new_bound_fluid_mass_fraction - in.composition[q][bound_fluid_idx];
                  double porosity_change = new_volume_fraction_porosity - volume_fraction_porosity;

                  // Do not allow negative bound water
                  if (bound_fluid_change + in.composition[q][bound_fluid_idx] < 0)
                    bound_fluid_change = -in.composition[q][bound_fluid_idx];
                  // Do not allow negative free water
                  if (porosity_change + volume_fraction_porosity < 0)
                    porosity_change = -volume_fraction_porosity;
                  // Apply the reactions
                  for (unsigned int c=0; c<in.composition[q].size(); ++c)
                    {
                      if (c == bound_fluid_idx && this->get_timestep_number() > 0)
                        // Apply the mass fraction change to the bound fluid content, which is a mass fraction
                        reaction_rate_out->reaction_rates[q][c] = bound_fluid_change / reaction_time_step_size;
                      else if (c == porosity_idx && this->get_timestep_number() > 0)
                        // Apply the volume fraction change to the porosity, which is a volume fraction
                        reaction_rate_out->reaction_rates[q][c] = porosity_change / reaction_time_step_size;
                      else
                        reaction_rate_out->reaction_rates[q][c] = 0.0;
                    }
                }
            }
        }
      else
        {
          katz2003_model.calculate_reaction_rate_outputs(in, out);
          katz2003_model.calculate_fluid_outputs(in, out, reference_T);
        }
    }

    template <int dim>
    void
    ReactiveFluidTransport<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Reactive Fluid Transport Model");
        {
          prm.enter_subsection("Katz 2003 model");
          {
            // read in melting model parameters
            ReactionModel::Katz2003MantleMelting<dim>::declare_parameters(prm);
          }
          prm.leave_subsection();

          prm.enter_subsection("Tian 2019 model");
          {
            ReactionModel::Tian2019Solubility<dim>::declare_parameters(prm);
          }
          prm.leave_subsection();

          prm.declare_entry("Base model","visco plastic",
                            Patterns::Selection(MaterialModel::get_valid_model_names_pattern<dim>()),
                            "The name of a material model incorporating the "
                            "addition of fluids. Valid values for this parameter "
                            "are the names of models that are also valid for the "
                            "``Material models/Model name'' parameter. See the documentation for "
                            "that for more information.");
          prm.declare_entry ("Reference fluid density", "2500",
                             Patterns::Double (0),
                             "Reference density of the melt/fluid$\\rho_{f,0}$. Units: \\si{\\kilogram\\per\\meter\\cubed}.");
          prm.declare_entry ("Shear to bulk viscosity ratio", "0.1",
                             Patterns::Double (0),
                             "Ratio between shear and bulk viscosity at the reference "
                             "permeability $\\phi_0=0.05$. The bulk viscosity additionally "
                             "scales with $\\phi_0/\\phi$. The shear viscosity is read in "
                             "from the base model. Units: dimensionless.");
          prm.declare_entry ("Minimum compaction viscosity", "0",
                             Patterns::Double (0),
                             "Lower cutoff for the compaction viscosity. Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Maximum compaction viscosity",
                             boost::lexical_cast<std::string>(std::numeric_limits<double>::max()),
                             Patterns::Double (0),
                             "Upper cutoff for the compaction viscosity. Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Reference fluid viscosity", "10",
                             Patterns::Double (0),
                             "The value of the constant melt/fluid viscosity $\\eta_f$. Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Exponential fluid weakening factor", "27",
                             Patterns::Double (0),
                             "The porosity dependence of the viscosity. Units: dimensionless.");
          prm.declare_entry ("Reference permeability", "1e-8",
                             Patterns::Double(),
                             "Reference permeability of the solid host rock."
                             "Units: \\si{\\meter\\squared}.");
          prm.declare_entry ("Fluid compressibility", "0.0",
                             Patterns::Double (0),
                             "The value of the compressibility of the fluid. "
                             "Units: \\si{\\per\\pascal}.");
          prm.declare_entry ("Fluid reaction time scale for operator splitting", "1e3",
                             Patterns::Double (0),
                             "In case the operator splitting scheme is used, the porosity field can not "
                             "be set to a new equilibrium fluid fraction instantly, but the model has to "
                             "provide a reaction time scale instead. This time scale defines how fast fluid "
                             "release and absorption happen, or more specifically, the parameter defines the "
                             "time after which the deviation of the porosity from the free fluid fraction "
                             "that would be in equilibrium with the solid will be reduced to a fraction of "
                             "$1/e$. So if the fluid reaction time scale is small compared "
                             "to the time step size, the reaction will be so fast that the porosity is very "
                             "close to this equilibrium value after reactions are computed. Conversely, "
                             "if the fluid reaction time scale is large compared to the time step size, almost no "
                             "fluid release and absorption will occur."
                             "\n\n"
                             "Also note that the fluid reaction time scale has to be larger than or equal to the reaction "
                             "time step used in the operator splitting scheme, otherwise reactions can not be "
                             "computed. If the model does not use operator splitting, this parameter is not used. "
                             "Units: yr or s, depending on the ``Use years "
                             "in output instead of seconds'' parameter.");
          prm.declare_entry ("Fluid-solid reaction scheme", "no reaction",
                             Patterns::Selection("no reaction|zero solubility|tian approximation|katz2003"),
                             "Select what type of scheme to use for reactions between fluid and solid phases. "
                             "The current available options are models where no reactions occur between "
                             "the two phases, or the solid phase is insoluble (zero solubility) and all "
                             "of the bound fluid is released into the fluid phase, tian approximation "
                             "use polynomials to describe hydration and dehydration reactions for four different "
                             "rock compositions as defined in Tian et al., 2019, or the Katz et. al. 2003 mantle "
                             "melting model. If the Katz 2003 melting model is used, its parameters are declared "
                             "in its own subsection.");
          prm.declare_entry ("Reference temperature", "293.",
                             Patterns::Double (0.),
                             "The reference temperature $T_0$ for the katz2003 reaction model. "
                             "The reference temperature is used in both the density and "
                             "viscosity formulas of this model. Units: \\si{\\kelvin}.");
        }
        prm.leave_subsection();

      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    ReactiveFluidTransport<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Reactive Fluid Transport Model");
        {
          AssertThrow( prm.get("Base model") != "reactive fluid transport",
                       ExcMessage("You may not use ``reactive fluid transport'' "
                                  "as the base model for the reactive fluid transport "
                                  "model itself.") );

          reference_rho_f                   = prm.get_double ("Reference fluid density");
          shear_to_bulk_viscosity_ratio     = prm.get_double ("Shear to bulk viscosity ratio");
          max_compaction_visc               = prm.get_double ("Maximum compaction viscosity");
          min_compaction_visc               = prm.get_double ("Minimum compaction viscosity");
          eta_f                             = prm.get_double ("Reference fluid viscosity");
          reference_permeability            = prm.get_double ("Reference permeability");
          alpha_phi                         = prm.get_double ("Exponential fluid weakening factor");
          fluid_compressibility             = prm.get_double ("Fluid compressibility");
          fluid_reaction_time_scale         = prm.get_double ("Fluid reaction time scale for operator splitting");
          reference_T                       = prm.get_double ("Reference temperature");

          // Create the base model and initialize its SimulatorAccess base
          // class; it will get a chance to read its parameters below after we
          // leave the current section.
          base_model = create_material_model<dim>(prm.get("Base model"));
          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(base_model.get()))
            sim->initialize_simulator (this->get_simulator());

          if (this->convert_output_to_years() == true)
            fluid_reaction_time_scale *= year_in_seconds;

          // Reaction scheme parameter
          if (prm.get ("Fluid-solid reaction scheme") == "zero solubility")
            {
              fluid_solid_reaction_scheme = zero_solubility;
            }
          else if (prm.get ("Fluid-solid reaction scheme") == "no reaction")
            {
              fluid_solid_reaction_scheme = no_reaction;
            }
          else if (prm.get ("Fluid-solid reaction scheme") == "tian approximation")
            {
              fluid_solid_reaction_scheme = tian_approximation;

              prm.enter_subsection("Tian 2019 model");
              {
                tian2019_model.initialize_simulator (this->get_simulator());
                tian2019_model.parse_parameters(prm);
              }
              prm.leave_subsection();
            }
          else if (prm.get ("Fluid-solid reaction scheme") == "katz2003")
            {
              fluid_solid_reaction_scheme = katz2003;
              prm.enter_subsection("Katz 2003 model");
              {
                katz2003_model.initialize_simulator (this->get_simulator());
                katz2003_model.parse_parameters(prm);
              }
              prm.leave_subsection();
            }
          else
            AssertThrow(false, ExcMessage("Not a valid fluid-solid reaction scheme"));

          if (fluid_solid_reaction_scheme == no_reaction)
            {
              AssertThrow(this->get_parameters().use_operator_splitting == false,
                          ExcMessage("The Fluid-reaction scheme no reaction should not be used with operator splitting."));
            }

          if (fluid_solid_reaction_scheme == zero_solubility)
            {
              AssertThrow(this->get_parameters().use_operator_splitting,
                          ExcMessage("The Fluid-reaction scheme zero solubility must be used with operator splitting."));
            }

          if (fluid_solid_reaction_scheme == tian_approximation)
            {
              AssertThrow(this->get_parameters().use_operator_splitting && this->get_parameters().reaction_solver_type == Parameters<dim>::ReactionSolverType::fixed_step,
                          ExcMessage("The Fluid-reaction scheme tian approximation must be used with operator splitting "
                                     "and 'Reaction solver type = fixed step'."));
            }

          if (this->get_parameters().use_operator_splitting)
            {
              if (this->get_parameters().reaction_solver_type == Parameters<dim>::ReactionSolverType::fixed_step)
                AssertThrow(fluid_reaction_time_scale >= this->get_parameters().reaction_time_step,
                            ExcMessage("The reaction time step " + Utilities::to_string(this->get_parameters().reaction_time_step)
                                       + " in the operator splitting scheme is too large to compute fluid release rates! "
                                       "You have to choose it in such a way that it is smaller than the 'Fluid reaction time scale for "
                                       "operator splitting' chosen in the material model, which is currently "
                                       + Utilities::to_string(fluid_reaction_time_scale) + "."));
              AssertThrow(fluid_reaction_time_scale > 0,
                          ExcMessage("The Fluid reaction time scale for operator splitting must be larger than 0!"));
            }

          AssertThrow(this->introspection().compositional_name_exists("porosity"),
                      ExcMessage("Material model Reactive Fluid Transport only "
                                 "works if there is a compositional field called porosity."));

          if (fluid_solid_reaction_scheme != katz2003)
            {
              AssertThrow(this->introspection().compositional_name_exists("bound_fluid"),
                          ExcMessage("Material model Reactive Fluid Transport only "
                                     "works if there is a compositional field called bound_fluid."));
            }
          else
            {
              AssertThrow(this->introspection().compositional_name_exists("peridotite"),
                          ExcMessage("Material model Katz 2003 Mantle Melting only "
                                     "works if there is a compositional field called peridotite."));
            }
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // After parsing the parameters for this model, parse parameters related to the base model.
      base_model->parse_parameters(prm);
      this->model_dependence = base_model->get_model_dependence();
      if (fluid_solid_reaction_scheme == zero_solubility)
        {
          AssertThrow(this->get_material_model().is_compressible() == false,
                      ExcMessage("The Fluid-reaction scheme zero solubility must be used with an incompressible base model."));
        }
    }



    template <int dim>
    void
    ReactiveFluidTransport<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (this->get_parameters().use_operator_splitting
          && out.template has_additional_output_object<ReactionRateOutputs<dim>>() == false)
        {
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::ReactionRateOutputs<dim>> (out.n_evaluation_points(), this->n_compositional_fields()));
        }
      base_model->create_additional_named_outputs(out);
    }
  }
}
// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(ReactiveFluidTransport,
                                   "reactive fluid transport",
                                   "Material model that is designed to advect fluids and compute "
                                   "fluid release and absorption based on different models for "
                                   "fluid-rock interaction. At present, models where no fluid-rock "
                                   "interactions occur or the solid has zero solubility are available. "
                                   "The properties of the solid can be taken from another material model "
                                   "that is used as a base model.")
  }
}
